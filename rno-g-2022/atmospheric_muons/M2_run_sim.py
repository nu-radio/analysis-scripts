import argparse
# import detector simulation modules
import NuRadioReco.modules.efieldToVoltageConverter
import NuRadioReco.modules.trigger.simpleThreshold
import NuRadioReco.modules.trigger.highLowThreshold
import NuRadioReco.modules.phasedarray.triggerSimulator
import NuRadioReco.modules.channelResampler
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.channelGenericNoiseAdder
import NuRadioReco.modules.RNO_G.hardwareResponseIncorporator
import NuRadioReco.modules.triggerTimeAdjuster
from NuRadioReco.utilities import units
import numpy as np
from NuRadioMC.simulation import simulation
import os
import copy

dipole_channels = [0]
dipole_noise_channels = [1]
deep_channels = [1, 2, 3, 4]
shallow_channels = [5, 6, 7]
channel_ids = [0, 1, 2, 3, 4, 5, 6, 7]
thresholds_pa = {'100Hz': 30.68, '1Hz': 38.62, '1mHz': 50.53}

passband_low = {}
passband_low_trigger = {}
passband_high = {}
filter_type = {}
order_low = {}
order_high = {}
passband_low_per_channel = {}
passband_high_per_channel = {}
for channel_id in deep_channels:
    passband_low_per_channel[channel_id] = [0, 220 * units.MHz]
    passband_high_per_channel[channel_id] = [96 * units.MHz, 100 * units.GHz]
    passband_low[channel_id] = [0 * units.MHz, 1000 * units.MHz]
    passband_low_trigger[channel_id] = [0 * units.MHz, 220 * units.MHz]
    passband_high[channel_id] = [96 * units.MHz, 100 * units.GHz]
    filter_type[channel_id] = 'cheby1'
    order_low[channel_id] = 7
    order_high[channel_id] = 4

for channel_id in dipole_channels:
    passband_low_per_channel[channel_id] = [0, 220 * units.MHz]
    passband_high_per_channel[channel_id] = [96 * units.MHz, 100 * units.GHz]
    passband_low[channel_id] = [0 * units.MHz, 1000 * units.MHz]
    passband_low_trigger[channel_id] = [0 * units.MHz, 220 * units.MHz]
    passband_high[channel_id] = [96 * units.MHz, 100 * units.GHz]
    filter_type[channel_id] = 'cheby1'
    order_low[channel_id] = 7
    order_high[channel_id] = 4

for channel_id in shallow_channels:
    passband_low_per_channel[channel_id] = [0 * units.MHz, 100 * units.GHz]
    passband_high_per_channel[channel_id] = [0 * units.MHz, 100 * units.GHz]
    passband_low[channel_id] = [1 * units.MHz, 1000 * units.MHz]
    passband_low_trigger[channel_id] = [1 * units.MHz, 180 * units.MHz]
    passband_high[channel_id] = [80 * units.MHz, 800 * units.GHz]
    filter_type[channel_id] = 'butter'
    order_low[channel_id] = 10
    order_high[channel_id] = 5

sampling_rate_pa = 0.5 * units.GHz
# x2 for upsampling
window_4ant = int(16 * units.ns * sampling_rate_pa * 2.0)
step_4ant = int(8 * units.ns * sampling_rate_pa * 2.0)
main_low_angle = np.deg2rad(-59.54968597864437)
main_high_angle = np.deg2rad(59.54968597864437)
phasing_angles_4ant = np.arcsin(np.linspace(np.sin(main_low_angle), np.sin(main_high_angle), 11))

efieldToVoltageConverter = NuRadioReco.modules.efieldToVoltageConverter.efieldToVoltageConverter()
simpleThreshold = NuRadioReco.modules.trigger.simpleThreshold.triggerSimulator()
highLowThreshold = NuRadioReco.modules.trigger.highLowThreshold.triggerSimulator()
phasedArrayTrigger = NuRadioReco.modules.phasedarray.triggerSimulator.triggerSimulator()
channelResampler = NuRadioReco.modules.channelResampler.channelResampler()
channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()
channelGenericNoiseAdder = NuRadioReco.modules.channelGenericNoiseAdder.channelGenericNoiseAdder()
hardwareResponseIncorporator = NuRadioReco.modules.RNO_G.hardwareResponseIncorporator.hardwareResponseIncorporator()

triggerTimeAdjuster = NuRadioReco.modules.triggerTimeAdjuster.triggerTimeAdjuster()
triggerTimeAdjuster.begin(pre_trigger_time=200 * units.ns)


class mySimulation(simulation.simulation):

    def _detector_simulation_filter_amp(self, evt, station, det):
        # hardware response
        hardwareResponseIncorporator.run(evt, station, det, sim_to_data=True)

    def _detector_simulation_trigger(self, evt, station, det):
        station_copy = copy.deepcopy(station)

        # 2) apply additional lowpass filter
        channelBandPassFilter.run(evt, station_copy, det,
                                  passband=passband_low_trigger, filter_type=filter_type, order=order_low, rp=0.1)

        # 3) apply trigger filter
        channelBandPassFilter.run(evt, station_copy, det, passband=passband_low_trigger, filter_type=filter_type, order=order_low, rp=0.1)
        channelBandPassFilter.run(evt, station_copy, det, passband=passband_high, filter_type=filter_type, order=order_high, rp=0.1)
        channelBandPassFilter.run(evt, station_copy, det, passband=passband_low, filter_type=filter_type, order=order_low, rp=0.1)

        Vrms_per_channel_copy = copy.deepcopy(self._Vrms_per_channel)
        print('Vrms_per_channel_copy[station_copy.get_id()]', Vrms_per_channel_copy[station_copy.get_id()], flush=True)

        # LPDAs
        threshold_high = {}
        threshold_low = {}
        for channel_id in shallow_channels:
            threshold_high[channel_id] = 2 * Vrms_per_channel_copy[station_copy.get_id()][channel_id]
            threshold_low[channel_id] = -2 * Vrms_per_channel_copy[station_copy.get_id()][channel_id]

        highLowThreshold.run(evt, station_copy, det,
                             threshold_high=threshold_high,
                             threshold_low=threshold_low,
                             coinc_window=50 * units.ns,
                             triggered_channels=shallow_channels,
                             number_concidences=2,  # 2/3 majority logic
                             trigger_name='hilo_2of3_2sigma')

        print('shallow trace', np.max(np.abs(station_copy.get_channel(shallow_channels[0]).get_trace())), flush=True)
        print('shallow trace threshold', threshold_high, flush=True)

        threshold_high = {}
        threshold_low = {}
        for channel_id in shallow_channels:
            threshold_high[channel_id] = 2.5 * Vrms_per_channel_copy[station_copy.get_id()][channel_id]
            threshold_low[channel_id] = -2.5 * Vrms_per_channel_copy[station_copy.get_id()][channel_id]

        highLowThreshold.run(evt, station_copy, det,
                             threshold_high=threshold_high,
                             threshold_low=threshold_low,
                             coinc_window=50 * units.ns,
                             triggered_channels=shallow_channels,
                             number_concidences=2,  # 2/3 majority logic
                             trigger_name='hilo_2of3_2.5sigma')

        threshold_high = {}
        threshold_low = {}
        for channel_id in shallow_channels:
            threshold_high[channel_id] = 3 * Vrms_per_channel_copy[station_copy.get_id()][channel_id]
            threshold_low[channel_id] = -3 * Vrms_per_channel_copy[station_copy.get_id()][channel_id]

        highLowThreshold.run(evt, station_copy, det,
                             threshold_high=threshold_high,
                             threshold_low=threshold_low,
                             coinc_window=50 * units.ns,
                             triggered_channels=shallow_channels,
                             number_concidences=2,  # 2/3 majority logic
                             trigger_name='hilo_2of3_3sigma')

        Vrms_deep = Vrms_per_channel_copy[station_copy.get_id()][deep_channels[0]]
        print("Vrms_deep", Vrms_deep, flush=True)

        trace = station_copy.get_channel(dipole_channels[0]).get_trace()
        print('dipole trace', np.max(np.abs(trace)), flush=True)
        print('dipole trace threshold', 1.5 * Vrms_deep, flush=True)

        # noiseless dipole
        simpleThreshold.run(evt, station_copy, det,
                            threshold=1.5 * Vrms_deep,
                            triggered_channels=dipole_channels,
                            trigger_name=f'dipole_1.5sigma')

        simpleThreshold.run(evt, station_copy, det,
                            threshold=2 * Vrms_deep,
                            triggered_channels=dipole_channels,
                            trigger_name=f'dipole_2sigma')

        simpleThreshold.run(evt, station_copy, det,
                            threshold=2.5 * Vrms_deep,
                            triggered_channels=dipole_channels,
                            trigger_name=f'dipole_2.5sigma')

        # noise dipole
        simpleThreshold.run(evt, station_copy, det,
                            threshold=2 * Vrms_deep,
                            triggered_channels=dipole_noise_channels,
                            trigger_name=f'deep_simple_2sigma')

        simpleThreshold.run(evt, station_copy, det,
                            threshold=2.5 * Vrms_deep,
                            triggered_channels=dipole_noise_channels,
                            trigger_name=f'deep_simple_2.5sigma')

        simpleThreshold.run(evt, station_copy, det,
                            threshold=3 * Vrms_deep,
                            triggered_channels=dipole_noise_channels,
                            trigger_name=f'deep_simple_3sigma')

        # PA
        phasedArrayTrigger.run(evt, station_copy, det,
                               Vrms=Vrms_deep,
                               threshold=thresholds_pa['100Hz'] * np.power(Vrms_deep, 2.0),
                               triggered_channels=deep_channels,
                               phasing_angles=phasing_angles_4ant,
                               ref_index=1.76,
                               trigger_name=f'PA_100Hz',  # the name of the trigger
                               trigger_adc=True,  # Don't have a seperate ADC for the trigger
                               adc_output=f'voltage',  # output in volts
                               trigger_filter=None,
                               upsampling_factor=2,
                               window=window_4ant,
                               step=step_4ant)

        phasedArrayTrigger.run(evt, station_copy, det,
                               Vrms=Vrms_deep,
                               threshold=thresholds_pa['1Hz'] * np.power(Vrms_deep, 2.0),
                               triggered_channels=deep_channels,
                               phasing_angles=phasing_angles_4ant,
                               ref_index=1.76,
                               trigger_name=f'PA_1Hz',  # the name of the trigger
                               trigger_adc=True,  # Don't have a seperate ADC for the trigger
                               adc_output=f'voltage',  # output in volts
                               trigger_filter=None,
                               upsampling_factor=2,
                               window=window_4ant,
                               step=step_4ant)

        phasedArrayTrigger.run(evt, station_copy, det,
                               Vrms=Vrms_deep,
                               threshold=thresholds_pa['1mHz'] * np.power(Vrms_deep, 2.0),
                               triggered_channels=deep_channels,
                               phasing_angles=phasing_angles_4ant,
                               ref_index=1.76,
                               trigger_name=f'PA_1mHz',  # the name of the trigger
                               trigger_adc=True,  # Don't have a seperate ADC for the trigger
                               adc_output=f'voltage',  # output in volts
                               trigger_filter=None,
                               upsampling_factor=2,
                               window=window_4ant,
                               step=step_4ant)

        # if(station_copy.has_triggered()):
        #     print("TRIGGERED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", flush=True)
        # 5) set trigger attributes of original station
        for trigger in station_copy.get_triggers().values():
            # print(trigger)
            station.set_trigger(trigger)

        # this module cuts the trace to the record length of the detector
        triggerTimeAdjuster.run(evt, station, det)

        if not 'trigger_names' in self._mout_attrs:
            # Might have some files which never ever trigger on a hybrid station and would consequently produce output trigger structures (nevents, 2) instead of (nevents, 6)
            # and only the LPDA trigger names. Merging and Veff calculation is then complicated.
            # This hack makes sure all triggers are written to the output hdf5. CAVEAT: make sure to adjust if trigger names above are changed!!!
            self._mout_attrs['trigger_names'] = ['hilo_2of3_2sigma', 'hilo_2of3_2.5sigma', 'hilo_2of3_3sigma',
                                                 'dipole_1.5sigma', 'dipole_2sigma', 'dipole_2.5sigma',
                                                 'deep_simple_2sigma', 'deep_simple_2.5sigma', 'deep_simple_3sigma',
                                                 'PA_100Hz', 'PA_1Hz', 'PA_1mHz']


results_folder = '/lustre/fs22/group/radio/lpyras/muon_sim/M2_results'
if not os.path.exists(results_folder):
    os.mkdir(results_folder)

parser = argparse.ArgumentParser(description='Run NuRadioMC simulation')
parser.add_argument('--inputfilename', type=str, default='/lustre/fs22/group/radio/lpyras/muon_sim/M1_simulation_input/M1_18.75_19.00_cos_theta_0.5_0.4_n_10_part_1.hdf5',
                    help='path to NuRadioMC input event list')
parser.add_argument('--detectordescription', type=str, default='/afs/ifh.de/group/radio/scratch/lpyras/muon_sim'
                                                               '/detector_4PA.json',
                    help='path to file containing the detector description')
parser.add_argument('--config', type=str, default='/afs/ifh.de/group/radio/scratch/lpyras/muon_sim/config.yaml',
                    help='NuRadioMC yaml config file')
parser.add_argument('--outputfilename', type=str, default=os.path.join(results_folder, 'M2_18.75_19.00_cos_theta_0.5_0.4_n_10_part_1.hdf5'),
                    help='hdf5 output filename')
parser.add_argument('--outputfilename_NuRadioReco', type=str, nargs='?', default=None,
                    help='outputfilename of NuRadioReco detector sim file')
args = parser.parse_args()


sim = mySimulation(inputfilename=args.inputfilename,
                   outputfilename=args.outputfilename,
                   detectorfile=args.detectordescription,
                   outputfilenameNuRadioReco=args.outputfilename_NuRadioReco,
                   config_file=args.config)
sim.run()
