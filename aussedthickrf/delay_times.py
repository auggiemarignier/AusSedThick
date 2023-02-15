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


def get_delay(stream: rf.RFStream()) -> rf.RFStream():
    """
    Gets the delay time between the first peak and the
    P-arrival (onset)
    """
    d = stream.times()[np.argmax(stream.data)] - (
        stream.statstream.onset - stream.stats.starttime
    )
    stream.stats["delay"] = d
    return stream


def main(rf_stream):
    rf_stream = rf_stream.select(channel="??R")

    # Create a dictionary of RFs by station
    rf_station_dict = defaultdict(list)
    for trc in rf_stream:
        rf_station_dict[trc.stats.station].append(trc)

    # Stack and compute delays
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
            s = get_delay(s)
            stacked.append(s)
            delays.append(s.stats["delay"])
    stacked = rf.RFStream(stacked).sort(["delay"], reverse=True)

    # Save delays
    stations_with_delays = {}
    with open(f"{basename}_delays.txt", "w") as f:
        for tr in stacked:
            stations_with_delays[tr.meta.station] = tr.stats.delay
            f.write(f"{tr.meta.station:<8}\t{tr.stats.delay:.2}\n")

    return stations_with_delays
