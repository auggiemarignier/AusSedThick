import requests
import pandas as pd
import geopandas as gpd
from os import path
import rf
import pygmt
import numpy as np
from scipy.signal import argrelmax, argrelmin


def get_twtt(autocorrelation: rf.rfstream.RFTrace) -> float:
    """
    Picks out the first local minimum of a radial RF autocorrelation.
    The autocorrelation should start at t=0 i.e. the latter half of the full autocorrelation
    """
    ind = argrelmin(autocorrelation.data)[0][0]
    return ind / autocorrelation.stats.sampling_rate


def get_tpsb(trace: rf.rfstream.RFTrace) -> float:
    """
    Finds the time of the first local maximum of trace.
    The pre-arrival part of the trace will have loads of little local maxima.
    Need to get the first maximum post P-onset
    """
    inds = argrelmax(trace.data)[0]
    p_arrival_ind = int(
        trace.stats.sampling_rate * (trace.stats.onset - trace.stats.starttime)
    )
    ind = inds[np.searchsorted(inds >= p_arrival_ind, True)]
    return ind / trace.stats.sampling_rate - (trace.stats.onset - trace.stats.starttime)


def get_geological_timeline() -> dict:
    """
    Extracts geological timeline data from the interactive geological timescale.

    Returns: (dict) geological_timeline
    """

    base_url = "https://raw.githubusercontent.com/i-c-stratigraphy/interactive-geological-timescale/master/src/assets/"

    r = requests.get(base_url + "timeline_data.json")
    _rjson = r.json()
    geological_timeline = {m["id"].split("/")[-1]: m for m in _rjson}

    r = requests.get(base_url + "time_interval_data.json")
    intervals = r.json()

    for K, V in geological_timeline.items():
        for k, v in intervals.items():
            if k.split("/")[-1] == K:
                V |= v
                break
    for K, V in geological_timeline.items():
        del V["name"]
        for k, v in V.items():
            if k == "narrow":
                V["children"] = [n.split("/")[-1] for n in v]
                del V[k]
            elif k == "broad":
                V["parent"] = [b.split("/")[-1] for b in v]
                del V[k]
                break

    return geological_timeline


def get_australian_sedimentary_basins() -> gpd.GeoDataFrame:
    """
    Extracts sedimentary basins from GA Portal
    Appends the period in a manner interoperable with the geological timeline

    Returns: (GeoDataBase) gdb
    """

    url = "https://services.ga.gov.au/gis/services/Australian_Geological_Provinces/MapServer/WFSServer"
    _eras = [
        "Cenozoic",
        "EarlyPaleozoic",
        "EarlyToLatePaleozoic",
        "LatePaleozoic",
        "Mesoproterozoic",
        "Mesozoic",
        "MesozoicCenozoic",
        "NeoproterozoicPaleozoic",
        "Paleoproterozoic",
        "PaleozoicCenozoic",
        "PaleozoicMesozoic",
    ]
    base_typeName = "Australian_Geological_Provinces:SedimentaryBasins-"
    params = {
        "request": "GetFeature",
        "service": "WFS",
        "version": "2.0.0",
        "typeName": base_typeName,
        "outputFormat": "geojson",
    }
    frames = []
    for era in _eras:
        p = path.join("..", "data", "australian_sedimentary_basins", era + ".json")
        try:
            d = gpd.read_file(p)
        except FileNotFoundError:
            params["typeName"] = base_typeName + era
            r = requests.get(url=url, params=params)
            if not r.ok:
                raise requests.exceptions.RequestException
            d = gpd.GeoDataFrame.from_features(r.json()["features"])
            d.to_file(p, driver="GeoJSON")
        except requests.exceptions.RequestException:
            raise LookupError(
                f"Not able to get sedimentary basin information for {era} era from file or GA Portal."
            )
        finally:
            frames.append(d)

    gdf = pd.concat(frames)
    gdf.set_index("GmlID", inplace=True)
    gdf.drop(gdf[gdf["onshoreOffshore"] == "Off"].index, inplace=True)

    geological_timeline = get_geological_timeline()
    periods = {}
    for k, v in geological_timeline.items():
        if v["type"] == "period":
            if v["parent"][0] in ["Cenozoic", "Mesozoic", "Paleozoic"]:
                periods[k] = v

    _periods = []
    fills = []
    ages = (
        []
    )  # mid point of period 'hasBeginning' and 'hasEnd', just to assign a numeric value
    for i, polygon in gdf.iterrows():
        current_len = len(fills)
        name = polygon["olderNameAge"]
        if not isinstance(name, str):
            print(i)
            continue
        for p in name.split():
            if p in geological_timeline:
                _p = p
                _type = geological_timeline[_p]["type"]
                going_up = True
                while _type != "period":
                    if going_up:
                        try:
                            _p = geological_timeline[_p]["parent"][0]
                        except IndexError:  # no more parents
                            going_up = False
                    else:
                        try:
                            _p = geological_timeline[_p]["children"][-1]  # oldest child
                        except IndexError as e:  # no more children
                            raise e(f"No period found for {p}")
                    _type = geological_timeline[_p]["type"]
                try:
                    fills.append(periods[_p]["fill"])
                    ages.append(
                        (periods[_p]["hasEnd"] + periods[_p]["hasBeginning"]) / 2
                    )
                    _periods.append(_p)
                except KeyError:  # precambrian
                    _periods.append("Precambrian")
                    fills.append("#80808080")
                    ages.append(2500)
                break
        if len(fills) == current_len:
            _periods.append(None)
            fills.append(None)
            ages.append(None)

    series = pd.DataFrame(
        {"fill": fills, "period_age": ages, "period": _periods}, gdf.index
    )
    gdf = pd.concat([gdf, series], axis=1)
    gdf.dropna(inplace=True, subset=["period_age"])
    return gdf


def australia_basemap(fig=None, frame=True, basins=True) -> pygmt.Figure:
    region = [112, 155, -46, -8]
    ln_min, ln_max, lt_min, lt_max = region
    projection = (
        f"M{int(np.mean([ln_min, ln_max]))}/{int(np.mean([lt_min, lt_max]))}/15c"
    )
    if fig is None:
        fig = pygmt.Figure()
    with pygmt.config(FONT_TITLE="20p,Helvetica,black"):
        fig.basemap(
            region=region,
            projection=projection,
            frame=frame,
        )
    fig.coast(
        region=region,
        projection=projection,
        shorelines=1,
        land="#ffffe6",
        water=None if basins else "#e6ffff",
    )
    if basins:
        gdf = get_australian_sedimentary_basins()
        gdf.set_index("provinceName", inplace=True)
        key_basins = gdf.loc[
            [
                "Eucla Basin",
                "Murray Basin",
                "Eromanga Basin",
                "Perth Basin",
                "McArthur Basin",
                "Amadeus Basin",
                "Southern Carnarvon Basin",
                "Kimberley Basin",
                "Canning Basin",
                "Karumba Basin",
                "Surat Basin",
                "Sydney Basin",
                "Gippsland Basin",
                "Otway Basin",
            ]
        ]
        key_basins = (
            key_basins.sort_values(by="period_age", ascending=False)
            .buffer(0.25)
            .buffer(-0.5)
            .buffer(0.25)
        )  # buffers smooth out the outlines
        fig.plot(data=key_basins, region=region, projection=projection, pen="1p,grey")
        fig.coast(  # replot water to hide offshore basins
            region=region,
            projection=projection,
            shorelines=1,
            resolution="i",
            water="#e6ffff",
        )
    return fig, region, projection