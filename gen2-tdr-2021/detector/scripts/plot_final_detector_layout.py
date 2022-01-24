from NuRadioReco.detector import generic_detector
from NuRadioReco.utilities import units
from matplotlib import pyplot as plt
import sys, os
import numpy as np

filename = sys.argv[1]
basename = os.path.splitext(os.path.basename(filename))[0]
det = generic_detector.GenericDetector(json_filename=filename)

xx = []
yy = []
station_ids = []
for station_id in det.get_station_ids():
    pos = det.get_absolute_position(station_id)
    xx.append(pos[0])
    yy.append(pos[1])
    station_ids.append(station_id)
station_ids = np.array(station_ids)
xx = np.array(xx)
yy = np.array(yy)

fig, ax = plt.subplots(1, 1)
ax.plot(xx, yy, "o")
ax.set_aspect('equal')
fig.savefig(f"{basename}.png")


print(basename)
dd = []
dd_deep = []
n_shallow = 0
n_hybrid = 0
mask_deep = station_ids < 2000
for i, station_id in enumerate(det.get_station_ids()):
    if(station_id < 2000):
        n_hybrid += 1
        dd_deep.append(np.sort(((xx[mask_deep] - xx[i])**2 + (yy[mask_deep] - yy[i])**2)**0.5)[1])
    else:
        n_shallow += 1
    dd.append(np.sort(((xx - xx[i])**2 + (yy - yy[i])**2)**0.5)[1])
dd_deep = np.unique(np.round(np.array(dd_deep)))
dd = np.unique(np.round(np.array(dd)))

print(f"number of stations: {len(xx)}")
print(f"number of hybrid stations: {n_hybrid}")
print(f"number of shallow stations: {n_shallow}")
print(f"station spacing shallow: {dd}")
print(f"station spacing deep: {dd_deep}")

p = np.array([xx, yy]).T

from shapely.geometry.polygon import Polygon
from scipy.spatial import Delaunay
tri = Delaunay(p)

coord_groups = [tri.points[x] for x in tri.simplices]
polygons = [Polygon(x) for x in coord_groups]

A = 0
for p in polygons:
    A += p.area
print(f"surface area = {A/units.km**2:.0f} km^2")