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



hybrid_center_spacing = 2.0
hyrbid_hexlength = 0.5
shallow_hexlength=0.75
nstations = 100

# define a basis vector
b0 = np.array([1,0])
b1 = np.array([np.cos(np.pi/3), np.sin(np.pi/3)])

b0scaled = b0 * hybrid_center_spacing * 1000. * 0.5
b1scaled = b1 * hybrid_center_spacing * 1000. * 0.5

# add points along the grid
intpoints = []
midpoints = []
size = 25.
for i in np.arange(-size, size+1,):
    for j in np.arange(-size, size+1):
        pt = b0 * i + b1 * (j - i)
        intpoints.append(pt)
        midpoints.append(pt + b0 * 0.5 - b1 * 0.5)
        midpoints.append(pt + b0 * 0.5)
        midpoints.append(pt + b1 * 0.5)

intpoints = np.array(intpoints) * hybrid_center_spacing * 1000.
midpoints = np.array(midpoints) * hybrid_center_spacing * 1000.
#plt.plot(intpoints[:,0], intpoints[:,1], 'bo')
#plt.plot(midpoints[:,0], midpoints[:,1], 'r.')

xx_deep_rot = intpoints.T[0]
yy_deep_rot = intpoints.T[1]
xx_shallow_rot = midpoints.T[0]
yy_shallow_rot = midpoints.T[1]

# cut out any parts of the array that don't fit in the dark sector wedge
max_rad = 20.5e3
min_h_rad = 2.5e3
max_h_rad = max_rad - 0.5 * hybrid_center_spacing*1e3

min_s_rad = 2.5e3
max_s_rad = max_rad

xx_deep_trim, yy_deep_trim = helper.softtrim_to_fit(xx_deep_rot, yy_deep_rot, min_h_rad, max_h_rad, b0scaled, b1scaled)
xx_shallow_trim, yy_shallow_trim = helper.trim_to_fit(xx_shallow_rot, yy_shallow_rot, min_s_rad, max_s_rad)

# adjust deep array
xx_deep_trim, yy_deep_trim = helper.adjust_stations(2.0*b0scaled, 2.0*b1scaled, xx_deep_trim, yy_deep_trim, [0,0], b0scaled, b1scaled)

#get shallow ring
xx_shallow_ring, yy_shallow_ring = helper.get_shallow_ring(b0scaled, b1scaled, xx_deep_trim, yy_deep_trim, xx_shallow_trim, yy_shallow_trim)
xx_shallow_ring_trim, yy_shallow_ring_trim = helper.softtrim_to_radius(xx_shallow_ring, yy_shallow_ring, min_s_rad, max_s_rad, b0scaled, b1scaled)


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
                "pos_altitude": None,
                "pos_site": "southpole",
                "station_id": int(number),
                "reference_station": 1001}

for it, number in enumerate(np.arange(2001, 2001+n_shallow_trim, 1)):
    dic[int(number)] = {
                "commission_time": "{TinyDate}:2017-11-04T00:00:00",
                "decommission_time": "{TinyDate}:2038-01-01T00:00:00",
                "pos_easting": int(xx_shallow_trim[it]-center_xx),
                "pos_northing": int(yy_shallow_trim[it]-center_yy),
                "pos_altitude": None,
                "pos_site": "southpole",
                "station_id": int(number),
                "reference_station": 2001}

for it, number in enumerate(np.arange(2001+n_shallow_trim, 2001+n_shallow_trim+n_shallow_ring_trim, 1)):
    dic[int(number)] = {
                "commission_time": "{TinyDate}:2017-11-04T00:00:00",
                "decommission_time": "{TinyDate}:2038-01-01T00:00:00",
                "pos_easting": int(xx_shallow_ring_trim[it]-center_xx),
                "pos_northing": int(yy_shallow_ring_trim[it]-center_yy),
                "pos_altitude": None,
                "pos_site": "southpole",
                "station_id": int(number),
                "reference_station": 2001}

with open('hex_shallow_array.json', 'w') as outfile:
    json.dump(dic, outfile, cls=NumpyEncoder, indent=4)