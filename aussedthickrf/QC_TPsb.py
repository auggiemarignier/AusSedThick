import os
import matplotlib.pyplot as plt
import tqdm
import numpy as np
from rf import read_rf, RFStream, IterMultipleComponents

from utils import calc_cc_conv, max_p_onset


time_delay = []
station = []
for station_5G in station_name:
    print("Doing", station_5G)
    try:
        data = os.path.join("rf_data/{}".format(station_5G), "")
        catfile = data + "rf_profile_events.xml"
        datafile = data + "rf_profile_data.h5"
        rffile = data + "rf_profile_rfs_1hz.h5"
        rffile_pdelay = data + "rf_profile_rfs_1Hz_pdel.h5"
        rffile_ref = data + "rf_profile_rfs_ref_1Hz_rvr.h5"
        #############################################################

        data_rf = read_rf(rffile, "H5")
        data = read_rf(datafile, "H5")
        print("No. of RF after SNR =", len(data_rf) / 3)

        data_rf_cc = RFStream()
        (cc_conv_all, data_rf_cc) = calc_cc_conv(data_rf, data)
        stream_rf = RFStream()
        for stream3c in tqdm(IterMultipleComponents(data_rf_cc, "onset", 3)):
            stream3c.taper(max_percentage=0.05)
            # stream3c.filter('bandpass', freqmin=0.1, freqmax=1)
            max_rf_r_abs = max(abs(stream3c.select(component="R")[0].data))
            max_rf_r = max(stream3c.select(component="R")[0].data)
            cc_conv = stream3c.select(component="R")[0].stats.cc_conv
            stream3c.select(component="R")[0].stats.max_rf_r_abs = max_rf_r_abs
            (max_abs_P, max_P) = max_p_onset(stream3c.select(component="R")[0])
            # if max of Direct P in R-rf is equal to the max in rest of signal
            if max_P == max_rf_r and max_P < 1:  # and cc_conv>0.4
                stream_rf.extend(stream3c)

        stream_rf.write(rffile_pdelay, "H5")
        print("No. of RF after direct P corrc. =", len(stream_rf) / 3)

        # Calculate Tpsb from Radial stack
        stack = stream_rf.stack()
        stack_R = stack.select(component="R")[0].data
        stack_Z = stack.select(component="Z")[0].data

        time_P = (np.argmax(stack_R) - np.argmax(stack_Z)) / stack[
            0
        ].stats.sampling_rate

        print("direct-P Time delay for station", station_5G, "=", time_P, "sec")
        time_delay.append(time_P)
        station.append(station_5G)

        kw = {
            "trim": (-2.5, 25),
            "fig_width": 6,
            "fillcolors": ("teal", "darkgrey"),
            "trace_height": 0.1,
            "show_vlines": "True",
            "scale": 2,
        }
        stream_rf.select(component="R", station="{}".format(station_5G)).sort(
            ["back_azimuth"]
        ).plot_rf(**kw)
        plt.savefig(
            "{}_{}_.1_1_directPD.pdf".format(station_5G, "R"),
            bbox_inches="tight",
            pad_inches=0.2,
        )
        # plt.show()
        plt.close()

    except:
        print("Data/rf not found for station", station_5G)
