import argparse
from NuRadioMC.simulation import simulation
import NuRadioReco.modules.efieldToVoltageConverter
import NuRadioReco.modules.trigger.simpleThreshold
import NuRadioReco.modules.phasedarray.triggerSimulator
import NuRadioReco.modules.channelBandPassFilter
from NuRadioReco.utilities import units
import numpy as np
import scipy
from scipy import constants
import matplotlib.pyplot as plt
import logging
import copy
import yaml

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger("runMB")

# 4 channel, 2x sampling, fft upsampling, 16 ns window
# 100 Hz -> 1.88
# 10 Hz -> 2.11
# 1 Hz -> 2.34

# 8 channel, 4x sampling, fft upsampling, 16 ns window
# 100 Hz -> 1.95
# 10 Hz -> 2.18
# 1 Hz -> 2.41

# initialize detector sim modules
simpleThreshold = NuRadioReco.modules.trigger.simpleThreshold.triggerSimulator()
channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()
phasedArrayTrigger = NuRadioReco.modules.phasedarray.triggerSimulator.triggerSimulator()

# assuming that PA consists out of 8 antennas (channel 0-7)
main_low_angle = np.deg2rad(-59.54968597864437)
main_high_angle = np.deg2rad(59.54968597864437)
phasing_angles_4ant = np.arcsin(np.linspace(np.sin(main_low_angle), np.sin(main_high_angle), 11))
phasing_angles_8ant = np.arcsin(np.linspace(np.sin(main_low_angle), np.sin(main_high_angle), 21))

passband_low = {}
passband_high = {}
filter_type = {}
order_low = {}
order_high = {}
for channel_id in range(0, 9):
    passband_low[channel_id] = [96 * units.MHz, 100 * units.GHz]
    passband_high[channel_id] = [0 * units.MHz, 220 * units.MHz]
    filter_type[channel_id] = 'cheby1'
    order_low[channel_id] = 4
    order_high[channel_id] = 7
passband_low[9] = [1 * units.MHz, 660 * units.MHz]
filter_type[9] = 'cheby1'
order_low[9] = 7
passband_high[9] = [100 * units.MHz, 100 * units.GHz]
order_high[9] = 4

class mySimulation(simulation.simulation):

    def _detector_simulation_filter_amp(self, evt, station, det):

        channelBandPassFilter.run(evt, station, det,
                                  passband=passband_low, filter_type=filter_type, order=order_low, rp=0.1)
        channelBandPassFilter.run(evt, station, det,
                                  passband=passband_high, filter_type=filter_type, order=order_high, rp=0.1)

    def _detector_simulation_trigger(self, evt, station, det):

        # channel 8 is a noiseless channel at 100m depth
        simpleThreshold.run(evt, station, det,
                                     threshold=3 * self._Vrms_per_channel[station.get_id()][8],
                                     triggered_channels=[8],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_3.0sigma')
        simpleThreshold.run(evt, station, det,
                                     threshold=2.5 * self._Vrms_per_channel[station.get_id()][8],
                                     triggered_channels=[8],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_2.5sigma')
        simpleThreshold.run(evt, station, det,
                                     threshold=2.0 * self._Vrms_per_channel[station.get_id()][8],
                                     triggered_channels=[8],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_2.0sigma')
        simpleThreshold.run(evt, station, det,
                                     threshold=1.5 * self._Vrms_per_channel[station.get_id()][8],
                                     triggered_channels=[8],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_1.5sigma')
        simpleThreshold.run(evt, station, det,
                                     threshold=1.0 * self._Vrms_per_channel[station.get_id()][8],
                                     triggered_channels=[8],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_1.0sigma')


        Vrms = self._Vrms_per_channel[station.get_id()][4]

        # run the 8 phased trigger
        # x4 for upsampling
        window_8ant = int(16 * units.ns * self._sampling_rate_detector * 4.0)
        step_8ant = int(8 * units.ns * self._sampling_rate_detector * 4.0)

        phasedArrayTrigger.run(evt, station, det,
                               Vrms = Vrms,
                               threshold = 1.95 * np.power(Vrms, 2.0) * window_8ant, # see phased trigger module for explanation                               
                               triggered_channels=range(0, 8),
                               phasing_angles=phasing_angles_8ant,
                               ref_index = 1.75,
                               trigger_name=f'PA_8channel_100Hz', # the name of the trigger
                               trigger_adc=False, # Don't have a seperate ADC for the trigger
                               adc_output=f'voltage', # output in volts
                               trigger_filter=None,
                               upsampling_factor=4,
                               window=window_8ant,
                               step=step_8ant)

        # run the 4 phased trigger
        # x2 for upsampling
        window_4ant = int(16 * units.ns * self._sampling_rate_detector * 2.0) 
        step_4ant = int(8 * units.ns * self._sampling_rate_detector * 2.0)

        phasedArrayTrigger.run(evt, station, det,
                               Vrms = Vrms,
                               threshold = 1.88 * np.power(Vrms, 2.0) * window_4ant,
                               triggered_channels=range(2, 6),
                               phasing_angles=phasing_angles_4ant,
                               ref_index = 1.75,
                               trigger_name=f'PA_4channel_100Hz', # the name of the trigger
                               trigger_adc=False, # Don't have a seperate ADC for the trigger
                               adc_output=f'voltage', # output in volts
                               trigger_filter=None,
                               upsampling_factor=2,
                               window=window_4ant,
                               step = step_4ant)

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
