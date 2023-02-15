from seismic import extract_event_traces
from seismic.receiver_fn import generate_rf, bulk_rf_report
import delay_times
from utils import signoise, parse_config

import rf
import sys
import yaml


try:
    config_file = sys.argv[1]
except IndexError:
    config_file = "config.yaml"

with open(config_file, "r") as f:
    config = yaml.safe_load(f)

#
# Data Collection
#

# Initial run to get EQ catalog
d = config["download_waveforms"]
extract_event_traces.main(
    catalog_only=True,
    inventory_file=d["inventory_file"],
    network_list=d["network_list"],
    station_list=d["station_list"],
    waveform_database=d["waveform_database"],
    event_catalog_file=d["event_catalog_file"],
    event_trace_datafile=d["event_trace_datafile"],
    start_time=d["start_time"],
    end_time=d["end_time"],
    taup_model=d["taup_model"],
    distance_range=d["distance_range"],
    magnitude_range=d["magnitude_range"],
    sw_magnitude_range=d["sw_magnitude_range"],
)

# Second run to get waveforms
extract_event_traces.main(
    catalog_only=False,
    inventory_file=d["inventory_file"],
    network_list=d["network_list"],
    station_list=d["station_list"],
    waveform_database=d["waveform_database"],
    event_catalog_file=d["event_catalog_file"],
    event_trace_datafile=d["event_trace_datafile"],
    start_time=d["start_time"],
    end_time=d["end_time"],
    taup_model=d["taup_model"],
    distance_range=d["distance_range"],
    magnitude_range=d["magnitude_range"],
    sw_magnitude_range=d["sw_magnitude_range"],
    resample_hz=d["resample_hz"],
    sw_resample_hz=d["sw_resample_hz"],
    p_data=d["p_data"],
    s_data=d["s_data"],
    sw_data=d["sw_data"],
)

#
# Calculate RFs
#
# some QC of waveforms is performed in generate_rf._main()
# can handle QC and preprocessing via config.json
generate_rf._main(
    input_file=None,
    output_file=None,
    network_list=None,
    station_list=None,
    config_file=None,
    only_corrections=None,
)

#
# RF QC
#
rf_stream_master = rf.read_rf(rf_filename)

# Drop RFs that do not meet quality estimate
MIN_SLOPE_RATIO = 5
rf_stream = rf.RFStream(
    [tr.copy() for tr in rf_stream_master if tr.stats.slope_ratio > MIN_SLOPE_RATIO]
).sort(["back_azimuth"])

rf_stream.taper(max_percentage=0.05)
rf_stream.trim2(-5, 10, reftime="onset")
rf_stream.moveout()


#
# Calcualte Delays
#
delay_times.main(rf_stream)


#
# Plots
#
bulk_rf_report.main()
