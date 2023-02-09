import sys
import os
import random
import numpy as np
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import rf
import rf.imaging

import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, FixedLocator, FixedFormatter, MaxNLocator)
import matplotlib
from matplotlib import cm
import pygmt
from collections import defaultdict


def plot_stacks(stream, fig_width=7., trace_height=0.5,
                stack_height=0.5, dpi=None,
                scale=1, trim=None,
                show_vlines=False):
    """
    Plot receiver functions.

    :param stream: stream to plot
    :param fname: filename to save plot to. Can be None. In this case
        the figure is left open.
    :param fig_width: width of figure in inches
    :param trace_height: height of one trace in inches
    :param stack_height: height of stack axes in inches
    :param dpi: dots per inch for the created figure
    :param scale: scale for individual traces
    :param fillcolors: fill colors for positive and negative wiggles
    :param trim: trim stream relative to onset before plotting using
         `~.rfstream.RFStream.slice2()`
    :param info: Plot one additional axes showing maximal two entries of
        the stats object. Each entry in this list is a list consisting of
        three entries: key, label and color.
        info can be None. In this case no additional axes is plotted.
    :param show_vlines: If True, show vertical alignment grid lines on plot
        at positions of the major x-tick marks.
    """

    if len(stream) == 0:
        return
    if trim:
        stream = stream.slice2(*trim, reftime='onset')
    N = len(stream)
    # calculate lag times
    stats = stream[0].stats
    times = stream[0].times() - (stats.onset - stats.starttime)
    # calculate axes and figure dimensions
    # big letters: inches, small letters: figure fraction
    H = trace_height
    HS = stack_height
    FB = 0.5
    FT = 0.2
    DW = 0.1
    FH = H * (N + 2) + HS + FB + FT + DW
    h = H / FH
    hs = HS / FH
    fb = FB / FH
    ft = FT / FH
    FL = 0.5
    FR = 0.2
    FW = fig_width
    FW3 = 0.8
    FW2 = FW - FL - FR
    fl = FL / FW
    fr = FR / FW
    fw2 = FW2 / FW
    fw3 = FW3 / FW
    # init figure and axes
    fig = plt.figure(figsize=(FW, FH), dpi=dpi)
    ax1 = fig.add_axes([fl, fb, fw2, h * (N + 2)])

    def _plot(ax, t, d, i, color, label):
        c1, c2 = (color, 'k')
        ax.text(-3, i + 0.2, label)
        if c1:
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor=c1)
        if c2:
            ax.fill_between(t, d + i, i, where=d < 0, lw=0., facecolor=c2)
        ax.plot(t, d + i, 'k')

    max_ = np.array([np.max(np.abs(tr.data)) for tr in stream])

    delays = np.array([tr.stats.delay for tr in stream])
    norm = matplotlib.colors.Normalize(vmin=np.min(delays), 
                                       vmax=np.max(delays))

    for i, tr in enumerate(stream):
        rgba_color = cm.turbo(tr.stats.delay)
        _plot(ax1, times, tr.data / max_[i] * scale, i + 1,
              rgba_color,
              f"{tr.stats.network}.{tr.stats.station}")

    # set x and y limits
    ax1.set_xlim(times[0], times[-1])
    ax1.set_ylim(-0.5, N + 1.5)
    ax1.set_yticklabels('')
    ax1.set_xlabel('time (s)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())

    sm = plt.cm.ScalarMappable(cmap='turbo', norm=norm)

    cbar = plt.colorbar(sm, orientation='horizontal')
    cbar.set_label('Delay [s]')
    return fig


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

client = Client("IRIS")
starttime = UTCDateTime("2001-01-01")
endtime = UTCDateTime("2023-01-02")
inventory = client.get_stations(
    network="AU,II,IU,G",
    starttime=starttime,
    endtime=endtime,
    level="channel",
    minlongitude=ln_min,
    maxlongitude=ln_max,
    minlatitude=lt_min,
    maxlatitude=lt_max,
)

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
mapfile = os.path.splitext(f"{basename}_delays_map.png")[0] + "_map.png"
fig.savefig(mapfile)
