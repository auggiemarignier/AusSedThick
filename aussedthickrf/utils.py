import requests
import pandas as pd
import geopandas as gpd


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
        params["typeName"] = base_typeName + era
        r = requests.get(url=url, params=params)
        frames.append(gpd.GeoDataFrame.from_features(r.json()["features"]))

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
    ages = []  # mid point of period 'hasBeginning' and 'hasEnd', just to assign a numeric value
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
                    fills.append(periods[_p]['fill'])
                    ages.append((periods[_p]['hasEnd'] + periods[_p]["hasBeginning"]) / 2)
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

    series = pd.DataFrame({"fill": fills, "period_age": ages, "period": _periods}, gdf.index)
    gdf = pd.concat([gdf, series], axis=1)
    gdf.dropna(inplace=True, subset=["period_age"])
    return gdf
