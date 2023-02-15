from seismic import extract_event_traces
from seismic.receiver_fn import generate_rf, bulk_rf_report, rf_quality_filter
import delay_times
from utils import signoise, parse_config

import rf
from os.path import splitext, expanduser
import sys
import json


try:
    config_file = sys.argv[1]
except IndexError:
    config_file = "config.json"

with open(config_file, "r") as f:
    config = json.load(f)

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
d = config["rfs"]
generate_rf._main(
    input_file=d["input_file"],
    output_file=d["output_file"],
    network_list=d["network_list"],
    station_list=d["station_list"],
    config_file=d["config_file"],
    only_corrections=d["only_corrections"],
)


#
# Quality filters
#
split = splitext(config["rfs"]["output_file"])
output_file = split[0] + "_qf" + split[1]
rf_quality_filter.main(
    input_file=config["rfs"]["output_file"],
    output_file=output_file,
    temp_dir=f"{expanduser('~')}/tmp",
    parallel=True,
)

#
# Calcualte Delays
#
delay_times.main(output_file)


#
# Plots
#
split = splitext(config["rfs"]["output_file"])
report_file = split[0] + "_bulkreport.pdf"
bulk_rf_report.main(
    output_file,
    report_file,
    network_list=config["download_waveforms"]["network_list"],
    station_list=config["download_waveforms"]["station_list"],
)
