from NuRadioReco.detector import generic_detector
from matplotlib import pyplot as plt
import sys, os

filename = sys.argv[1]
basename = os.path.splitext(os.path.basename(filename))[0]
print(basename)
det = generic_detector.GenericDetector(json_filename=filename)

xx = []
yy = []
for station_id in det.get_station_ids():
    pos = det.get_absolute_position(station_id)
    xx.append(pos[0])
    yy.append(pos[1])

fig, ax = plt.subplots(1, 1)
ax.plot(xx, yy, "o")
ax.set_aspect('equal')
fig.savefig(f"{basename}.png")