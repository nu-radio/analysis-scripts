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
import NuRadioReco.modules.phasedarray.triggerSimulator
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
phasedArrayTrigger = NuRadioReco.modules.phasedarray.triggerSimulator.triggerSimulator()

triggerTimeAdjuster = NuRadioReco.modules.triggerTimeAdjuster.triggerTimeAdjuster()

thresholds = {
  '2/4_100Hz': 3.9498194908011524,
  '2/4_10mHz': 4.919151494949084,
  '2/6_100Hz': 4.04625348733533,
  '2/6_10mHz': 5.015585491483261,
  'fhigh': 0.15,
  'flow': 0.08
  }

main_low_angle = np.deg2rad(-59.54968597864437)
main_high_angle = np.deg2rad(59.54968597864437)
phasing_angles_4ant = np.arcsin(np.linspace(np.sin(main_low_angle), np.sin(main_high_angle), 11))
phasing_angles_8ant = np.arcsin(np.linspace(np.sin(main_low_angle), np.sin(main_high_angle), 21))

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
for channel_id in range(4, 13):
    passband_low[channel_id] = [96 * units.MHz, 100 * units.GHz]
    passband_high[channel_id] = [0 * units.MHz, 220 * units.MHz]
    filter_type[channel_id] = 'cheby1'
    order_low[channel_id] = 4
    order_high[channel_id] = 7


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

        # channel 8 is a noiseless channel at 100m depth
        simpleThreshold.run(evt, station, det,
                                     threshold=4 * self._Vrms_per_channel[station.get_id()][12],
                                     triggered_channels=[12],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_4sigma')
        simpleThreshold.run(evt, station, det,
                                     threshold=3 * self._Vrms_per_channel[station.get_id()][12],
                                     triggered_channels=[12],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_3.0sigma')
        simpleThreshold.run(evt, station, det,
                                     threshold=2.5 * self._Vrms_per_channel[station.get_id()][12],
                                     triggered_channels=[12],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_2.5sigma')
        simpleThreshold.run(evt, station, det,
                                     threshold=2.0 * self._Vrms_per_channel[station.get_id()][12],
                                     triggered_channels=[12],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_2.0sigma')
        simpleThreshold.run(evt, station, det,
                                     threshold=1.5 * self._Vrms_per_channel[station.get_id()][12],
                                     triggered_channels=[12],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_1.5sigma')
        simpleThreshold.run(evt, station, det,
                                     threshold=1.0 * self._Vrms_per_channel[station.get_id()][12],
                                     triggered_channels=[12],  # run trigger on all channels
                                     number_concidences=1,
                                     trigger_name=f'dipole_1.0sigma')

        Vrms = self._Vrms_per_channel[station.get_id()][8]
        # run the 8 phased trigger
        # x4 for upsampling
        window_8ant = int(16 * units.ns * self._sampling_rate_detector * 4.0)
        step_8ant = int(8 * units.ns * self._sampling_rate_detector * 4.0)

        phasedArrayTrigger.run(evt, station, det,
                               Vrms=Vrms,
                               threshold=62.15 * np.power(Vrms, 2.0),  # see phased trigger module for explanation
                               triggered_channels=range(4, 12),
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
                               threshold=30.85 * np.power(Vrms, 2.0),
                               triggered_channels=range(6, 10),
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

sim = mySimulation(inputfilename=args.inputfilename,
                            outputfilename=args.outputfilename,
                            detectorfile=args.detectordescription,
                            outputfilenameNuRadioReco=args.outputfilenameNuRadioReco,
                            config_file=args.config,
                            default_detector_station=1)
sim.run()

