import os
from tqdm import tqdm
import matplotlib.pyplot as plt
from obspy import read_events
from obspy.clients.fdsn import Client
from rf import read_rf, RFStream, iter_event_data, IterMultipleComponents

from utils import signoise


client = Client("IRIS")
client_5g = Client("http://auspass.edu.au:80")  # or for any temporary network
inventory = client_5g.get_stations(
    network="7I",
    station="*",
    level="response",
    minlatitude=-39,
    maxlatitude=-25,
    minlongitude=128,
    maxlongitude=141.5,
)

station_names = [
    f"{network.code}.{station.code}" for network in inventory for station in network
]

# loop for getting Rf's station by station

for station_name in station_names:
    print("Doing", station_name)
    network, station = station_name.split(".")
    datapath = os.path.join("rf_data", station_name, "")
    catfile = datapath + "rf_profile_events.xml"
    datafile = datapath + "rf_profile_data.h5"
    rffile = datapath + "rf_profile_rfs_1hz.h5"
    plotpath = "plots"

    if not os.path.exists(datapath):
        os.makedirs(datapath)
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)

    st_lat = inventory.get_coordinates(f"{network}.{station}..BHZ")[
        "latitude"
    ]  # change here for diff networks
    st_long = inventory.get_coordinates(f"{network}.{station}..BHZ")[
        "longitude"
    ]  # change here for diff networks
    station_coord = inventory.get_coordinates(f"{network}.{station}..BHZ")
    # add a try loop if stations have different data formate for example EH*, HH* etc.
    starttime = inventory.select(station=station)[0][0].start_date
    endtime = inventory.select(station=station)[0][0].end_date
    # next bit acquires the teleseismic earthquakes for a station from IRIS.
    if not os.path.exists(catfile):
         catalog = client.get_events(
            starttime=starttime,
            endtime=endtime,
            minmagnitude=5.5,
            latitude=st_lat,
            longitude=st_long,
            minradius=30,
            maxradius=95,
        )
        catalog.write(catfile, "QUAKEML")

    catalog = read_events(catfile)

    print("# of events:", len(catalog))

    if not os.path.exists(datafile):
        stream = RFStream()

        with tqdm() as pbar:
            # try:
            for s in iter_event_data(
                catalog,
                inventory.select(station=station),
                client_5g.get_waveforms,
                pbar=pbar,
            ):  # ,**kw):
                SNR = signoise(
                    s.select(component="Z")[0]
                )  # SNR is computed on Z component
                if SNR >= 1.5:
                    s.detrend("linear")
                    s.detrend("demean").taper(max_percentage=0.05)
                    stream.extend(s)

        stream.write(datafile, "H5")
        print(f"Len of data per component after SNR: {len(stream) // 3}")

    try:
        data = read_rf(datafile, "H5")
        stream = RFStream()

        for stream3c in tqdm(IterMultipleComponents(data, "onset", 3)):
            stream3c.detrend("linear")
            stream3c.detrend("demean").taper(max_percentage=0.05)
            stream3c.filter(
                "bandpass", freqmin=0.1, freqmax=1
            )  # Change frequency content here
            stream3c.trim2(-25, 75, "onset")
            if len(stream3c) != 3:
                continue
            # kw = {'gauss':.75}
            stream3c.rf(
                rotate="NE->RT", deconvolve="time"
            )  # change here for different deconvolution options.
            # stream3c.rf(rotate='NE->RT',deconvolve='iterative',**kw)

            stream3c.moveout()
            stream.extend(stream3c)

        stream.write(rffile, "H5")

        # plots obtained Rf's
        component = "R"
        kw = {
            "trim": (-2.5, 25),
            "fig_width": 6,
            "fillcolors": ("navy", "darkgrey"),
            "trace_height": 0.1,
            "show_vlines": "True",
            "scale": 3.5,
        }
        stream.select(component=component, station=station).sort(
            ["back_azimuth"]
        ).plot_rf(**kw)
        plt.savefig(
            f"{plotpath}/{station_name}_{component}_.1_1.pdf",
            bbox_inches="tight",
            pad_inches=0.2,
        )
        plt.close()

        print("No. of RF=", len(stream) // 3, "for station", station_name)
        print("--------------------------------------------------\n")
    except:
        print(f"No data for station {station_name}")
