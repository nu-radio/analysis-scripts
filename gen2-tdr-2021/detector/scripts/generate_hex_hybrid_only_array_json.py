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


hybrid_spacing = 1.75
nstations = 1000

## make a big rectangular grid of stations, one for deep and one for shallow
## these will extend outside the dark sector by a lot. we will trim them down later.
# xx_deep = np.arange(-44, 44.1, hybrid_spacing)*1000.
# yy_deep = np.arange(-44, 44.1, hybrid_spacing)*1000.
# d = np.meshgrid(xx_deep, yy_deep)
# xx_deep = d[0].flatten()
# yy_deep = d[1].flatten()
#
## rotate the points by four degrees (because it aligns them a bit better with axis somehow)
# xx_deep_rot = []
# yy_deep_rot = []
# xx_shallow_rot = []
# yy_shallow_rot = []
#
#
# for x, y in zip(xx_deep, yy_deep):
#     x1, y1 = helper.rotate_point(x,y,0,0,4)
#     xx_deep_rot.append(x1)
#     yy_deep_rot.append(y1)

hybrid_center_spacing = 1.75
hyrbid_hexlength = 0.5
shallow_hexlength=0.75
nstations = 100

# define a basis vector
b0 = np.array([1,0])
b1 = np.array([np.cos(np.pi/3), np.sin(np.pi/3)])

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


# cut out any parts of the array that don't fit in the dark sector wedge
min_rad = 2.5e3
max_rad = 23e3
xx_deep_trim, yy_deep_trim = helper.trim_to_fit(xx_deep_rot, yy_deep_rot, min_rad, max_rad)



center_xx = np.average(xx_deep_trim)
center_yy = np.average(yy_deep_trim)
print(center_xx, center_yy)

n_deep_trim = len(xx_deep_trim)


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

with open('hex_hybrid_only_array.json', 'w') as outfile:
    json.dump(dic, outfile, cls=NumpyEncoder, indent=4)