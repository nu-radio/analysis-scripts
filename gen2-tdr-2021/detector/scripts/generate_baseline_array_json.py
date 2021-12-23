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


# make up a geometry baseline array
hybrid_spacing = 1.75
shallow_spacing = 1.75
nstations = 1000

# define a basis vector
b0 = np.array([1,0])
b1 = np.array([1,1])
b0scaled = b0 * shallow_spacing * 1000.
b1scaled = b1 * shallow_spacing * 1000.
offset = 0.5*b1scaled

# make a big rectangular grid of stations, one for deep and one for shallow
# these will extend outside the dark sector by a lot. we will trim them down later.
xx_deep = np.arange(-44, 44.1, hybrid_spacing)*1000.
yy_deep = np.arange(-44, 44.1, hybrid_spacing)*1000.
d = np.meshgrid(xx_deep, yy_deep)
xx_deep = d[0].flatten()
yy_deep = d[1].flatten()

xx_shallow = np.arange(-44 - 0.5 * shallow_spacing, 44.1 + 0.5 * shallow_spacing, shallow_spacing)*1000.
yy_shallow = np.arange(-44 - 0.5 * shallow_spacing, 44.1 + 0.5 * shallow_spacing, shallow_spacing)*1000.
s = np.meshgrid(xx_shallow, yy_shallow)
xx_shallow = s[0].flatten()
yy_shallow = s[1].flatten()


# rotate the points by four degrees (because it aligns them a bit better with axis somehow)
xx_deep_rot = []
yy_deep_rot = []
xx_shallow_rot = []
yy_shallow_rot = []
for x, y in zip(xx_shallow, yy_shallow):
    x1, y1 = helper.rotate_point(x,y,0,0,4)
    xx_shallow_rot.append(x1)
    yy_shallow_rot.append(y1)

for x, y in zip(xx_deep, yy_deep):
    x1, y1 = helper.rotate_point(x,y,0,0,4)
    xx_deep_rot.append(x1)
    yy_deep_rot.append(y1)

# rotate scaled basis vectors
x1, y1 = helper.rotate_point(b0scaled[0], b0scaled[1],0,0,4)
b0scaled = np.array([x1,y1])
x1, y1 = helper.rotate_point(b1scaled[0], b1scaled[1],0,0,4)
b1scaled = np.array([x1,y1])

# cut out any parts of the array that don't fit in the dark sector wedge
min_rad = 2.5e3
max_rad = 22.e3
xx_deep_trim, yy_deep_trim = helper.trim_to_fit(xx_deep_rot, yy_deep_rot, min_rad, max_rad)
xx_shallow_trim, yy_shallow_trim = helper.trim_to_fit(xx_shallow_rot, yy_shallow_rot, min_rad, max_rad)

# adjust deep array
xx_deep_trim, yy_deep_trim = helper.adjust_stations(b0scaled, b1scaled, xx_deep_trim, yy_deep_trim, offset)
# get shallow ring
xx_shallow_ring, yy_shallow_ring = helper.get_shallow_ring(b0scaled, b1scaled, xx_deep_trim, yy_deep_trim, xx_shallow_trim, yy_shallow_trim, offset)
xx_shallow_ring_trim, yy_shallow_ring_trim = helper.softtrim_to_radius(xx_shallow_ring, yy_shallow_ring, min_rad, max_rad, b0scaled, b1scaled)


center_xx = np.average(np.concatenate([xx_deep_trim, xx_shallow_trim, xx_shallow_ring_trim]))
center_yy = np.average(np.concatenate([yy_deep_trim, yy_shallow_trim, yy_shallow_ring_trim]))
print(center_xx, center_yy)

n_deep_trim = len(xx_deep_trim)
n_shallow_trim = len(xx_shallow_trim)
n_shallow_ring_trim = len(xx_shallow_ring_trim)

dic = {}
for it, number in enumerate(np.arange(1001, 1001+n_deep_trim, 1)):
    dic[int(number)] = {
                "commission_time": "{TinyDate}:2017-11-04T00:00:00",
                "decommission_time": "{TinyDate}:2038-01-01T00:00:00",
                "pos_easting": int(xx_deep_trim[it]-center_xx),
                "pos_northing": int(yy_deep_trim[it]-center_yy),
                "pos_altitude":2800,
                "pos_site": "southpole",
                "station_id": int(number),
                "reference_station": 1001}

for it, number in enumerate(np.arange(2001, 2001+n_shallow_trim, 1)):
    dic[int(number)] = {
                "commission_time": "{TinyDate}:2017-11-04T00:00:00",
                "decommission_time": "{TinyDate}:2038-01-01T00:00:00",
                "pos_easting": int(xx_shallow_trim[it]-center_xx),
                "pos_northing": int(yy_shallow_trim[it]-center_yy),
                "pos_altitude":2800,
                "pos_site": "southpole",
                "station_id": int(number),
                "reference_station": 2001}

for it, number in enumerate(np.arange(2001+n_shallow_trim, 2001+n_shallow_trim+n_shallow_ring_trim, 1)):
    dic[int(number)] = {
                "commission_time": "{TinyDate}:2017-11-04T00:00:00",
                "decommission_time": "{TinyDate}:2038-01-01T00:00:00",
                "pos_easting": int(xx_shallow_ring_trim[it]-center_xx),
                "pos_northing": int(yy_shallow_ring_trim[it]-center_yy),
                "pos_altitude":2800,
                "pos_site": "southpole",
                "station_id": int(number),
                "reference_station": 2001}

with open('baseline_array.json', 'w') as outfile:
    json.dump(dic, outfile, cls=NumpyEncoder, indent=4)