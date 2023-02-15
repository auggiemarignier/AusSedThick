from seismic import extract_event_traces
from seismic.receiver_fn import generate_rf, bulk_rf_report
import delay_times
from utils import signoise

import rf

#
# Data Collection
#

# Initial run to get EQ catalog
extract_event_traces.main(
    inventory_file=None,
    network_list="*",
    station_list=None,
    waveform_database=None,
    event_catalog_file=None,
    event_trace_datafile=None,
    start_time=None,
    end_time=None,
    taup_model=None,
    distance_range=None,
    magnitude_range=None,
    sw_magnitude_range=None,
    catalog_only=True,
    resample_hz=None,
    sw_resample_hz=None,
    p_data=None,
    s_data=None,
    sw_data=None,
    dry_run=None,
)

# Second run to get waveforms
extract_event_traces.main(
    inventory_file=None,
    network_list=None,
    station_list=None,
    waveform_database=None,
    event_catalog_file=None,
    event_trace_datafile=None,
    start_time=None,
    end_time=None,
    taup_model=None,
    distance_range=None,
    magnitude_range=None,
    sw_magnitude_range=None,
    catalog_only=False,
    resample_hz=None,
    sw_resample_hz=None,
    p_data=None,
    s_data=None,
    sw_data=None,
    dry_run=None,
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
