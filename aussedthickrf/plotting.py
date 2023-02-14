
# Plot stacks
fig = plot_stacks(stacked)
fig.savefig(f"{basename}_delays.png")


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
