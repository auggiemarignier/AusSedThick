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

# Load precomputed RFs for a network
RF_FILENAME = sys.argv[1]
rf_stream_master = rf.read_rf(RF_FILENAME)
basename = os.path.splitext(RF_FILENAME)[0]

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

# Plot stacks
fig = plot_stacks(stacked)
fig.savefig(f"{basename}_delays.png")

# Save delays
stations_with_delays = {}
with open(f"{basename}_delays.txt", "w") as f:
    for tr in stacked:
        stations_with_delays[tr.meta.station] = tr.stats.delay
        f.write(f"{tr.meta.station:<8}\t{tr.stats.delay:.2}\n")


# Plot map
ln_min, ln_max = (112, 155)
lt_min, lt_max = (-46, -8)

try:
    inventory = read_inventory(sys.argv[2])
except IndexError:
    starttime = UTCDateTime("1950-01-01")
    endtime = UTCDateTime("2023-01-02")
    client = Client("IRIS")
    perm_inventory = client.get_stations(
        network="AU,II,IU,G",
        starttime=starttime,
        endtime=endtime,
        level="channel",
        minlongitude=ln_min,
        maxlongitude=ln_max,
        minlatitude=lt_min,
        maxlatitude=lt_max,
    )
    client = Client("AUSPASS")
    temp_inventory = client.get_stations(
        starttime=starttime,
        endtime=endtime,
        level="channel",
        minlongitude=ln_min,
        maxlongitude=ln_max,
        minlatitude=lt_min,
        maxlatitude=lt_max,
    )
    inventory = perm_inventory + temp_inventory

all_stations = list(
    {
        sta.code
        for network in inventory
        for sta in network
        if sta.code in stations_with_delays
    }  # using a set to avoid duplicates
)
lats = np.zeros_like(all_stations, dtype=float)
lons = np.zeros_like(all_stations, dtype=float)
nets = np.zeros_like(all_stations, dtype=str)
delays = np.zeros_like(all_stations, dtype=float)
for i, sta in enumerate(all_stations):
    network = inventory.select(station=sta)[0]
    station = network[0]
    lats[i] = station.latitude
    lons[i] = station.longitude
    nets[i] = network.code
    delays[i] = stations_with_delays[sta]

fig = pygmt.Figure()
fig.basemap(region=[ln_min, ln_max, lt_min, lt_max], frame=True)
fig.coast(shorelines=1, land="#ffffe6", water="#e6ffff", borders="2/1p,grey")

markers = "dhist"
pygmt.makecpt(cmap="turbo", series=[delays.min(), delays.max()])

marker = random.choice(markers)
fig.plot(
    x=lons,
    y=lats,
    style=f"{marker}c",
    fill=delays,
    cmap=True,
    size=np.full_like(lons, 0.5),
)
fig.colorbar(frame="af+lDelay Time TPsb (s)")
mapfile = f"{basename}_map.png"
fig.savefig(mapfile)
