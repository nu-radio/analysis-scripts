THRES=20
# ntries:
# a reasonable quantity is probably simulating 1 minute of detector livetime
# ntries = 60 sec * 2.4 GHz / 2048 samples ~= 7e7
python PREP01_hpol_PA_threshold_optimisation.py --detectordescription ../detector/single_station/HPol_PA_channels.json --ntries 70000000 --threshold $THRES
