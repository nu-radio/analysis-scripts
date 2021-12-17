import argparse
from NuRadioMC.simulation import simulation
import NuRadioReco.modules.efieldToVoltageConverter
import NuRadioReco.modules.trigger.simpleThreshold
import NuRadioReco.modules.phasedarray.triggerSimulator
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.efieldToVoltageConverter
import NuRadioReco.modules.channelAddCableDelay
import NuRadioReco.modules.channelResampler
import NuRadioReco.modules.channelGenericNoiseAdder
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


# 4 channel PA: 2x sampling, fft upsampling, 16 ns window, 11 beams 
# Linear fit: -0.2519 x + 0.727
# Quadratic fit: -0.0002581 x x - 0.2393 x + 0.5845
# Noise rate -> threshold
#
# 10000.000 Hz -> 22.74
# 100.000 Hz -> 30.68
# 1.000 Hz -> 38.62
# 0.001 Hz -> 50.53

# 8 Channel PA: 4x sampling, fft upsampling, 16 ns window, 21 beams
# Linear fit: -0.134 x + 1.295
# Quadratic fit: -0.0007308 x x - 0.06398 x - 0.2857
# Noise rate -> threshold
#
# 10000.000 Hz -> 46.98
# 100.000 Hz -> 61.90
# 1.000 Hz -> 76.83
# 0.001 Hz -> 99.22

# initialize detector sim modules
simpleThreshold = NuRadioReco.modules.trigger.simpleThreshold.triggerSimulator()
channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()
phasedArrayTrigger = NuRadioReco.modules.phasedarray.triggerSimulator.triggerSimulator()
channelAddCableDelay = NuRadioReco.modules.channelAddCableDelay.channelAddCableDelay()
efieldToVoltageConverter = NuRadioReco.modules.efieldToVoltageConverter.efieldToVoltageConverter()
channelResampler = NuRadioReco.modules.channelResampler.channelResampler()
channelGenericNoiseAdder = NuRadioReco.modules.channelGenericNoiseAdder.channelGenericNoiseAdder()

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


class TDR_Simulation(simulation.simulation):
    # TODO:
    # - ensure bandpass is applied on the fly and saved traces are without the passband
    # - might need to implement _part1 to apply passband on the fly to all channels?
    # - ensure Vrms is calculated correctly (caveat: passband)
    def _detector_simulation_part2(self):
        # start detector simulation
        efieldToVoltageConverter.run(self._evt, self._station, self._det)  # convolve efield with antenna pattern
        # downsample trace to internal simulation sampling rate (the efieldToVoltageConverter upsamples the trace to
        # 20 GHz by default to achive a good time resolution when the two signals from the two signal paths are added)
        channelResampler.run(self._evt, self._station, self._det, sampling_rate=1. / self._dt)
  
        if self._is_simulate_noise():
            max_freq = 0.5 / self._dt
            channel_ids = self._det.get_channel_ids(self._station.get_id())
            Vrms = {}
            for channel_id in channel_ids:
                norm = self._bandwidth_per_channel[self._station.get_id()][channel_id]
                Vrms[channel_id] = self._Vrms_per_channel[self._station.get_id()][channel_id] / (norm / (max_freq)) ** 0.5  # normalize noise level to the bandwidth its generated for
                channelGenericNoiseAdder.run(self._evt, self._station, self._det, amplitude=Vrms, min_freq=0 * units.MHz,
                                        max_freq=max_freq, type='rayleigh', excluded_channels=self._noiseless_channels[self._station])
  
        self._detector_simulation_filter_amp(self._evt, self._station, self._det)
        self._detector_simulation_trigger(self._evt, self._station, self._det)


    def _detector_simulation_filter_amp(self, evt, station, det):
        # in the filter amp we don't do anything, apply passband in the _detector_simulation_trigger instead (after calc. Vrms)
        pass

    def _detector_simulation_trigger(self, evt, station, det):

        # TODO do we need a noiseless channel at 150m depth? I think we don't
        #simple_thresholds = [] # list of n-sigma thresholds to run # TODO do not run any for now
        #for n_sigma in simple_thresholds:
        #    simpleThreshold.run(evt, station, det,
        #                             threshold=n_sigma * self._Vrms_per_channel[station.get_id()][8],
        #                             triggered_channels=[9],  # run trigger on defined channel # TODO: replace by a noiseless channel
        #                             number_concidences=1,
        #                             trigger_name=f'dipole_{n_sigma}sigma')

        # get the Vrms before applying the passband
        Vrms = self._Vrms_per_channel[station.get_id()][9]

        # apply the BandPassFilter just before the trigger
        channelBandPassFilter.run(evt, station, det,
                                passband=passband_low, filter_type=filter_type, order=order_low, rp=0.1)
        channelBandPassFilter.run(evt, station, det,
                                passband=passband_high, filter_type=filter_type, order=order_high, rp=0.1)        

        # run the 8 phased trigger
        # x4 for upsampling
        window_8ant = int(16 * units.ns * self._sampling_rate_detector * 4.0)
        step_8ant = int(8 * units.ns * self._sampling_rate_detector * 4.0)

        phasedArrayTrigger.run(evt, station, det,
                               Vrms=Vrms,
                               threshold=61.90 * np.power(Vrms, 2.0),  # see phased trigger module for explanation
                               triggered_channels=range(9, 9+8),
                               phasing_angles=phasing_angles_8ant,
                               ref_index=1.75,
                               trigger_name=f'PA_8channel_100Hz',  # the name of the trigger
                               trigger_adc=False,  # Don't have a seperate ADC for the trigger
                               adc_output=f'voltage',  # output in volts
                               trigger_filter=None,
                               upsampling_factor=4,
                               window=window_8ant,
                               step=step_8ant)

        # run the 4 phased trigger
        # x2 for upsampling
        window_4ant = int(16 * units.ns * self._sampling_rate_detector * 2.0)
        step_4ant = int(8 * units.ns * self._sampling_rate_detector * 2.0)

        phasedArrayTrigger.run(evt, station, det,
                               Vrms=Vrms,
                               threshold=30.68 * np.power(Vrms, 2.0),
                               triggered_channels=range(9, 9+4),
                               phasing_angles=phasing_angles_4ant,
                               ref_index=1.75,
                               trigger_name=f'PA_4channel_100Hz',  # the name of the trigger
                               trigger_adc=False,  # Don't have a seperate ADC for the trigger
                               adc_output=f'voltage',  # output in volts
                               trigger_filter=None,
                               upsampling_factor=2,
                               window=window_4ant,
                               step=step_4ant)


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

sim = TDR_Simulation(inputfilename=args.inputfilename,
                            outputfilename=args.outputfilename,
                            detectorfile=args.detectordescription,
                            outputfilenameNuRadioReco=args.outputfilenameNuRadioReco,
                            config_file=args.config,
                            default_detector_station=1)
sim.run()
