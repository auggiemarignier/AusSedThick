import numpy as np

ff = open("P_delay_lat_long.txt", "w")
for station_5G in station_name:
    st_lat = inventory.get_coordinates("7I.{}..BHZ".format(station_5G))[
        "latitude"
    ]  # change here for diff networks
    st_long = inventory.get_coordinates("7I.{}..BHZ".format(station_5G))[
        "longitude"
    ]  # change here for diff networks

    for i in range(len(time_delay)):
        if station_5G == station[i]:
            if time_delay[i] < 0.58:
                basement_depth = 366 * time_delay[i]  # for shallow sediments
            else:
                basement_depth = (
                    3206.9 * time_delay[i] - 1661.2
                )  # for thicker sediments

            ff.write(
                "{} {} {} {} {}\n".format(
                    st_lat,
                    st_long,
                    np.around(time_delay[i], 2),
                    basement_depth,
                    station[i],
                )
            )

ff.close()
