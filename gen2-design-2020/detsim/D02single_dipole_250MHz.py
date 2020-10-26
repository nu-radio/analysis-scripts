import argparse
# import detector simulation modules
import NuRadioReco.modules.trigger.highLowThreshold
import NuRadioReco.modules.trigger.simpleThreshold
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.triggerTimeAdjuster
from NuRadioReco.utilities import units
import yaml
import numpy as np
from scipy import constants
from NuRadioMC.simulation import simulation
import scipy
import matplotlib.pyplot as plt
import logging
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger("runMB")

# initialize detector sim modules
simpleThreshold = NuRadioReco.modules.trigger.simpleThreshold.triggerSimulator()
highLowThreshold = NuRadioReco.modules.trigger.highLowThreshold.triggerSimulator()
channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()

triggerTimeAdjuster = NuRadioReco.modules.triggerTimeAdjuster.triggerTimeAdjuster()


class mySimulation(simulation.simulation):

    def _detector_simulation_filter_amp(self, evt, station, det):

        channelBandPassFilter.run(evt, station, det,
                                  passband=[1 * units.MHz, 230 * units.MHz], filter_type="cheby1", order=9, rp=0.1)
        channelBandPassFilter.run(evt, station, det,
                                  passband=[80 * units.MHz, 800 * units.GHz], filter_type="cheby1", order=4, rp=0.1)

    def _detector_simulation_trigger(self, evt, station, det):
        threshold = {}
        for channel_id in det.get_channel_ids(station.get_id()):
            threshold[channel_id] = 2.5 * self._Vrms_per_channel[station.get_id()][channel_id]
        simpleThreshold.run(evt, station, det,
                                     threshold=threshold,
                                     triggered_channels=None,  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_2.5sigma')  # the name of the trigger
        threshold = {}
        for channel_id in det.get_channel_ids(station.get_id()):
            threshold[channel_id] = 1.5 * self._Vrms_per_channel[station.get_id()][channel_id]
        simpleThreshold.run(evt, station, det,
                                     threshold=threshold,
                                     triggered_channels=None,  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_1.5sigma')  # the name of the trigger
        threshold = {}
        for channel_id in det.get_channel_ids(station.get_id()):
            threshold[channel_id] = 1.0 * self._Vrms_per_channel[station.get_id()][channel_id]
        simpleThreshold.run(evt, station, det,
                                     threshold=threshold,
                                     triggered_channels=None,  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_1.0sigma')  # the name of the trigger
#         triggerTimeAdjuster.run(evt, station, det)


parser = argparse.ArgumentParser(description='Run NuRadioMC simulation')
parser.add_argument('inputfilename', type=str,
                    help='path to NuRadioMC input event list')
parser.add_argument('detectordescription', type=str,
                    help='path to file containing the detector description')
parser.add_argument('config', type=str,
                    help='NuRadioMC yaml config file')
parser.add_argument('outputfilename', type=str,
                    help='hdf5 output filename')
parser.add_argument('outputfilenameNuRadioReco', type=str, nargs='?', default=None,
                    help='outputfilename of NuRadioReco detector sim file')
args = parser.parse_args()

sim = mySimulation(inputfilename=args.inputfilename,
                            outputfilename=args.outputfilename,
                            detectorfile=args.detectordescription,
                            outputfilenameNuRadioReco=args.outputfilenameNuRadioReco,
                            config_file=args.config,
                            default_detector_station=1)
sim.run()

