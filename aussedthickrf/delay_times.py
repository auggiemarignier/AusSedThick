import numpy as np
import rf
from os.path import splitext
from collections import defaultdict


def rf_quality_control(stream: rf.RFStream) -> rf.RFStream():
    # Drop RFs that do not meet quality estimate
    MIN_SLOPE_RATIO = 5
    stream = rf.RFStream(
        [tr.copy() for tr in stream if tr.stats.slope_ratio > MIN_SLOPE_RATIO]
    ).sort(["back_azimuth"])

    stream.taper(max_percentage=0.05)
    stream.trim2(-5, 10, reftime="onset")
    stream.moveout()
    return stream


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


def main(input_file: str) -> rf.RFStream():
    rf_stream_master = rf.read_rf(input_file)

    rf_stream = rf_quality_control(rf_stream_master)
    rf_stream = rf_stream.select(channel="??R")

    # Create a dictionary of RFs by station
    rf_station_dict = defaultdict(list)
    for trc in rf_stream:
        rf_station_dict[trc.stats.station].append(trc)

    # Stack and compute delays
    stacked = []
    delays = []
    for k, v in rf_station_dict.items():
        try:
            # sometimes this fails
            # moveout might produce some weird numpy arrays that can't be stacked
            s = v.stack()[0]
        except ValueError:
            continue
        else:
            s = get_delay(s)
            stacked.append(s)
            delays.append(s.stats["delay"])
    stacked = rf.RFStream(stacked).sort(["delay"], reverse=True)

    # Save delays
    output_file = f"{splitext(input_file)[0]}_delays.txt"
    with open(output_file, "w") as f:
        for tr in stacked:
            f.write(f"{tr.meta.station:<8}\t{tr.stats.delay}\n")

    return stacked
