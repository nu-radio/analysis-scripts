import numpy as np
import itertools
import json
import matplotlib.pyplot as plt
import helper as helper

class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)



shallow_center_spacing = 1.0
hybrid_center_spacing = 3.0
nstations = 100

# define a basis vector
b0 = np.array([1,0])
b1 = np.array([np.cos(np.pi/3), np.sin(np.pi/3)])
b2 = np.array([0,1])
b3 = np.array([-np.sin(np.pi/3), np.cos(np.pi/3)])

b0scaled = b0 * shallow_center_spacing * 1000.
b1scaled = b1 * shallow_center_spacing * 1000.

# add points along the grid
intpoints = []
deeppoints = []
size = 50.
for i in np.arange(-size, size+1,):
    for j in np.arange(-size, size+1):
        pt = b0 * i + b1 * (j - i)
        intpoints.append(pt)
        deeppoints.append(pt)

intpoints = np.array(intpoints) * shallow_center_spacing * 1000.
deeppoints = np.array(deeppoints) * hybrid_center_spacing * 1000.
#plt.plot(intpoints[:,0], intpoints[:,1], 'bo')
#plt.plot(midpoints[:,0], midpoints[:,1], 'r.')

xx_deep_rot = deeppoints.T[0]
yy_deep_rot = deeppoints.T[1]
xx_shallow_rot = intpoints.T[0]
yy_shallow_rot = intpoints.T[1]

# cut out any parts of the array that don't fit in the dark sector wedge
max_rad = 23.9e3
min_h_rad = 2.5e3
max_h_rad = max_rad

min_s_rad = 2.5e3
max_s_rad = max_rad

xx_deep_trim, yy_deep_trim = helper.trim_to_fit(xx_deep_rot, yy_deep_rot, min_h_rad, max_h_rad)
xx_shallow_trim, yy_shallow_trim = helper.trim_to_fit(xx_shallow_rot, yy_shallow_rot, min_s_rad, max_s_rad)
xx_shallow_ring, yy_shallow_ring = helper.get_shallow_ring(b0scaled, b1scaled, xx_deep_trim, yy_deep_trim, xx_shallow_trim, yy_shallow_trim)
xx_shallow_ring_trim, yy_shallow_ring_trim = helper.softtrim_to_radius(xx_shallow_ring, yy_shallow_ring, min_s_rad, max_s_rad, b0scaled, b1scaled)


# filter out duplicates of deep and shallow positions
xx_shallow_trim_filtered = []
yy_shallow_trim_filtered = []
nclose = 0
for x, y in zip(xx_shallow_trim, yy_shallow_trim):
    if np.any(np.isclose(x, np.array(xx_deep_trim), atol=1e-2) & np.isclose(y, np.array(yy_deep_trim))):
        nclose += 1
        continue
    else:
        xx_shallow_trim_filtered.append(x)
        yy_shallow_trim_filtered.append(y)
print(f"skipped {nclose} shallow stations, because there are deep already at the same position")

center_xx = np.average(np.concatenate([xx_deep_trim, xx_shallow_trim_filtered, xx_shallow_ring_trim]))
center_yy = np.average(np.concatenate([yy_deep_trim, yy_shallow_trim_filtered, yy_shallow_ring_trim]))
print(center_xx, center_yy)

n_deep_trim = len(xx_deep_trim)
n_shallow_trim_filtered = len(xx_shallow_trim_filtered)
n_shallow_ring_trim = len(xx_shallow_ring_trim)

dic = {}
for it, number in enumerate(np.arange(1001, 1001+n_deep_trim, 1)):
    dic[int(number)] = {
                "commission_time": "{TinyDate}:2017-11-04T00:00:00",
                "decommission_time": "{TinyDate}:2038-01-01T00:00:00",
                "pos_easting": int(xx_deep_trim[it]-center_xx),
                "pos_northing": int(yy_deep_trim[it]-center_yy),
                "pos_altitude": None,
                "pos_site": "southpole",
                "station_id": int(number),
                "reference_station": 1001}

for it, number in enumerate(np.arange(2001, 2001+n_shallow_trim_filtered, 1)):
    dic[int(number)] = {
                "commission_time": "{TinyDate}:2017-11-04T00:00:00",
                "decommission_time": "{TinyDate}:2038-01-01T00:00:00",
                "pos_easting": int(xx_shallow_trim_filtered[it]-center_xx),
                "pos_northing": int(yy_shallow_trim_filtered[it]-center_yy),
                "pos_altitude": None,
                "pos_site": "southpole",
                "station_id": int(number),
                "reference_station": 2001}

for it, number in enumerate(np.arange(2001+n_shallow_trim_filtered, 2001+n_shallow_trim_filtered+n_shallow_ring_trim, 1)):
    dic[int(number)] = {
                "commission_time": "{TinyDate}:2017-11-04T00:00:00",
                "decommission_time": "{TinyDate}:2038-01-01T00:00:00",
                "pos_easting": int(xx_shallow_ring_trim[it]-center_xx),
                "pos_northing": int(yy_shallow_ring_trim[it]-center_yy),
                "pos_altitude": None,
                "pos_site": "southpole",
                "station_id": int(number),
                "reference_station": 2001}

with open('hex_shallowheavy_array.json', 'w') as outfile:
    json.dump(dic, outfile, cls=NumpyEncoder, indent=4)
