from NuRadioReco.detector import generic_detector
from matplotlib import pyplot as plt
import sys

det = generic_detector.GenericDetector(json_filename=sys.argv[1])

xx = []
yy = []
for station_id in det.get_station_ids():
    pos = det.get_absolute_position(station_id)
    xx.append(pos[0])
    yy.append(pos[1])

fig, ax = plt.subplots(1, 1)
ax.plot(xx, yy)
plt.show()