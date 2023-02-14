import sys
import os
import random
import numpy as np
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.core.inventory.inventory import read_inventory
import rf
import rf.imaging
import pygmt
from collections import defaultdict

from utils import plot_stacks


def main(rf_filename):
    # Load precomputed RFs for a network
    rf_stream_master = rf.read_rf(rf_filename)
    basename = os.path.splitext(rf_filename)[0]

    # Drop RFs that do not meet quality estimate
    MIN_SLOPE_RATIO = 5
    rf_stream = rf.RFStream([tr.copy() for tr in rf_stream_master
                            if tr.stats.slope_ratio > MIN_SLOPE_RATIO]).sort(['back_azimuth'])
    rf_stream = rf_stream.select(channel='??R')

    # Create a dictionary of RFs by station
    rf_station_dict = defaultdict(list)
    for trc in rf_stream:
        rf_station_dict[trc.stats.station].append(trc)

    # Trim RFs
    for k, v in rf_station_dict.items():
        temp_stream = rf.RFStream([tr.copy() for tr in v])
        rf_station_dict[k] = temp_stream.trim2(-5, 10, reftime='onset')

    # Moveout, stack and compute delays
    stacked = []
    delays = []
    for k, v in rf_station_dict.items():
        m = v.copy().moveout()
        try:
            # sometimes this fails
            # moveout might produce some weird numpy arrays that can't be stacked
            s = m.stack()[0]
        except ValueError:
            continue
        else:
            d = s.times()[np.argmax(s.data)] - (s.stats.onset - s.stats.starttime)
            s.stats['delay'] = d
            stacked.append(s)
            delays.append(d)
    stacked = rf.RFStream(stacked).sort(['delay'], reverse=True)

    # Save delays
    stations_with_delays = {}
    with open(f"{basename}_delays.txt", "w") as f:
        for tr in stacked:
            stations_with_delays[tr.meta.station] = tr.stats.delay
            f.write(f"{tr.meta.station:<8}\t{tr.stats.delay:.2}\n")

    return stations_with_delays

