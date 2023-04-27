import numpy as np
import matplotlib
from matplotlib import cm
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import requests
import pandas as pd
import geopandas as gpd


def signoise(self):  # , winsig, winnoise, relative='onset'):
    """
    Determine signal noise ratio by dividing the maximum in the two windows.
    """
    st = self.stats
    self_copy = self.copy()
    self_copy.detrend().taper(max_percentage=0.05)
    self_copy.filter("bandpass", freqmin=0.1, freqmax=1)  # ,corners=2, zerophase=True)
    winsig = [-5, 25]  # signal window
    winnoise = [-45, -15]  # noise window
    rel_time = getattr(st, "onset")
    winsig0 = rel_time - st.starttime + winsig[0]
    winsig1 = rel_time - st.starttime + winsig[1]
    winnoise0 = rel_time - st.starttime + winnoise[0]
    winnoise1 = rel_time - st.starttime + winnoise[1]
    #
    t = np.arange(self.stats.npts) * 1.0 / st.sampling_rate
    datasig = self_copy.data[(t >= winsig0) * (t <= winsig1)]
    datanoise = self_copy.data[(t >= winnoise0) * (t <= winnoise1)]

    try:
        st.signoise = max(abs(datasig)) / max(abs(datanoise))
        return st.signoise
    except:
        st.signoise = 0
        return st.signoise


def calc_cc_conv(data_rf, data):
    """
    Returns rf stream with correlation coeffecient in tr.stats.cc_conv
    """
    rf_all = data_rf
    sig_all = data.copy()
    rf_all.sort(keys=["onset"])
    sig_all.sort(keys=["onset"])
    sig_all.rotate("NE->RT")
    sig_all.trim2(-25, 75, "onset")
    ###
    cc = []
    fit = []
    if len(rf_all) == len(sig_all):
        # print('True')
        for i in range(0, int(len(rf_all) / 3)):
            # for i in range(0,2)):
            try:
                rf = rf_all[3 * i : 3 * i + 3]
                sig = sig_all[3 * i : 3 * i + 3]
                # rf,sig are stream containing same three traces for RF and ZNE
                obs_Z = sig.select(component="Z")[0].copy()
                obs_R = sig.select(component="R")[0].copy()
                obs_rfR = rf.select(component="R")[0].copy()
                sr = obs_Z.stats.sampling_rate
                # Filter using SNR bandpass
                obs_Z.detrend().taper(max_percentage=0.05, max_length=2.0)
                obs_R.detrend().taper(max_percentage=0.05, max_length=2.0)
                obs_Z.filter("bandpass", freqmin=0.1, freqmax=1.0, zerophase=True)
                obs_R.filter("bandpass", freqmin=0.1, freqmax=1.0, zerophase=True)
                obs_rfR.filter("bandpass", freqmin=0.1, freqmax=1.0, zerophase=True)

                pred_R = obs_R.copy()
                pred_R.stats.channel = "PRR"
                st = pred_R.stats.starttime
                time_shift = pred_R.stats.onset - pred_R.stats.starttime  # rel P onset
                ind1 = int(np.ceil(time_shift * sr))

                ind2 = ind1 + len(obs_Z.data)  # [leadin:leadin + n]
                pred_R.data = np.convolve(obs_Z.data, obs_rfR.data, mode="full")[
                    ind1:ind2
                ]

                obs_Z.trim(st + 20, st + 45)
                obs_R.trim(st + 20, st + 45)
                pred_R.trim(st + 20, st + 45)
                rf.select(component="R")[0].stats.cc_conv = np.corrcoef(
                    obs_R.data * obs_R.data, pred_R.data * pred_R.data
                )[0][1]
                # rf.select(component='R')[0].stats.cc_conv = np.corrcoef(abs(obs_R.data), abs(pred_R.data))[0][1]
                cc.append(np.corrcoef(obs_R.data, pred_R.data)[0][1])
            except:
                print("  ")  # do nothing
    return cc, rf_all


def max_p_onset(self):
    """
    Determine max amplitude of r-rf around P-onset. winP_st and winP_end should be adjusted in case of very thick sediments i.e. > 2km.
    """
    onset = self.stats.onset
    self_copy = self.copy()
    self_copy.trim(onset - 5, onset + 30)
    self_copy.taper(max_percentage=0.05)
    st = self_copy.stats
    rel_time = getattr(st, "onset")

    winP_st = rel_time - st.starttime - 0.1  # P onset -1
    winP_end = rel_time - st.starttime + 1.5  # P onset +1
    t = np.arange(st.npts) * 1.0 / st.sampling_rate
    max_abs = max(
        abs(self_copy.data[(t >= winP_st) * (t <= winP_end)])
    )  # gets max amp around P (-.1,1.5 sec)
    max_ = max(self_copy.data[(t >= winP_st) * (t <= winP_end)])
    self.stats.max_abs_P = max_abs
    self.stats.max_P = max_
    return max_abs, max_


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
