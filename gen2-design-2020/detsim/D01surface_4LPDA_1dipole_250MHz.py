import argparse
# import detector simulation modules
import NuRadioReco.modules.efieldToVoltageConverter
import NuRadioReco.modules.trigger.highLowThreshold
import NuRadioReco.modules.trigger.simpleThreshold
import NuRadioReco.modules.channelResampler
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.channelGenericNoiseAdder
import NuRadioReco.modules.triggerTimeAdjuster
import NuRadioReco.modules.channelAntennaDedispersion
import NuRadioReco.modules.ARIANNA.hardwareResponseIncorporator
import NuRadioReco.modules.custom.deltaT.calculateAmplitudePerRaySolution
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

thresholds = {
  '2/4_100Hz': 3.9498194908011524,
  '2/4_10mHz': 4.919151494949084,
  '2/6_100Hz': 4.04625348733533,
  '2/6_10mHz': 5.015585491483261,
  'fhigh': 0.15,
  'flow': 0.08
  }

passband_low = {}
passband_high = {}
filter_type = {}
order_low = {}
order_high = {}
for channel_id in range(0, 4):
    passband_low[channel_id] = [1 * units.MHz, thresholds['fhigh']]
    passband_high[channel_id] = [thresholds['flow'], 800 * units.GHz]
    filter_type[channel_id] = 'butter'
    order_low[channel_id] = 10
    order_high[channel_id] = 5
passband_low[4] = [1 * units.MHz, 230 * units.MHz]
filter_type[4] = 'cheby1'
order_low[4] = 9
passband_high[4] = [80 * units.MHz, 100 * units.GHz]
order_high[4] = 4


class mySimulation(simulation.simulation):

    def _detector_simulation_filter_amp(self, evt, station, det):
        channelBandPassFilter.run(evt, station, det,
                                  passband=passband_low, filter_type=filter_type, order=order_low, rp=0.1)
        channelBandPassFilter.run(evt, station, det,
                                  passband=passband_high, filter_type=filter_type, order=order_high, rp=0.1)

    def _detector_simulation_trigger(self, evt, station, det):
        # run a high/low trigger on the 4 downward pointing LPDAs
        threshold_high = {}
        threshold_low = {}
        for channel_id in det.get_channel_ids(station.get_id()):
            threshold_high[channel_id] = 2 * self._Vrms_per_channel[station.get_id()][channel_id]
            threshold_low[channel_id] = -2 * self._Vrms_per_channel[station.get_id()][channel_id]
        highLowThreshold.run(evt, station, det,
                                    threshold_high=threshold_high,
                                    threshold_low=threshold_low,
                                    coinc_window=40 * units.ns,
                                    triggered_channels=[0, 1, 2, 3],  # select the LPDA channels
                                    number_concidences=2,  # 2/4 majority logic
                                    trigger_name='LPDA_2of4_2sigma')

        threshold_high = {}
        threshold_low = {}
        for channel_id in det.get_channel_ids(station.get_id()):
            threshold_high[channel_id] = thresholds['2/4_100Hz'] * self._Vrms_per_channel[station.get_id()][channel_id]
            threshold_low[channel_id] = -thresholds['2/4_100Hz'] * self._Vrms_per_channel[station.get_id()][channel_id]
        highLowThreshold.run(evt, station, det,
                                    threshold_high=threshold_high,
                                    threshold_low=threshold_low,
                                    coinc_window=40 * units.ns,
                                    triggered_channels=[0, 1, 2, 3],  # select the LPDA channels
                                    number_concidences=2,  # 2/4 majority logic
                                    trigger_name='LPDA_2of4_100Hz',
                                    set_not_triggered=not self._station.has_triggered())

        threshold_high = {}
        threshold_low = {}
        for channel_id in det.get_channel_ids(station.get_id()):
            threshold_high[channel_id] = thresholds['2/4_10mHz'] * self._Vrms_per_channel[station.get_id()][channel_id]
            threshold_low[channel_id] = -thresholds['2/4_10mHz'] * self._Vrms_per_channel[station.get_id()][channel_id]
        highLowThreshold.run(evt, station, det,
                                    threshold_high=threshold_high,
                                    threshold_low=threshold_low,
                                    coinc_window=40 * units.ns,
                                    triggered_channels=[0, 1, 2, 3],  # select the LPDA channels
                                    number_concidences=2,  # 2/4 majority logic
                                    trigger_name='LPDA_2of4_10mHz',
                                    set_not_triggered=not self._station.has_triggered())
        simpleThreshold.run(evt, station, det,
                                     threshold=2.5 * self._Vrms_per_channel[station.get_id()][4],
                                     triggered_channels=[4],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name='dipole_2.5sigma')  # the name of the trigger
        simpleThreshold.run(evt, station, det,
                                     threshold=1.5 * self._Vrms_per_channel[station.get_id()][4],
                                     triggered_channels=[4],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name='dipole_1.5sigma')  # the name of the trigger
        simpleThreshold.run(evt, station, det,
                                     threshold=1 * self._Vrms_per_channel[station.get_id()][4],
                                     triggered_channels=[4],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name='dipole_1.0sigma')  # the name of the trigger
        triggerTimeAdjuster.run(evt, station, det)


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

