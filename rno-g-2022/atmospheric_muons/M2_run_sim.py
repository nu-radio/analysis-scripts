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
channel_ids = [0, 1, 2, 3, 4, 5, 6]
thresholds_pa = {}
thresholds_pa['1Hz'] = 38.62
thresholds_pa['1mHz'] = 50.53

passband_low_per_channel = {}
passband_high_per_channel = {}
for channel_id in deep_channels:
    passband_low_per_channel[channel_id] = {}
    passband_high_per_channel[channel_id] = {}
    passband_low_per_channel[channel_id] = [0, 220 * units.MHz]
    passband_high_per_channel[channel_id] = [96 * units.MHz, 100 * units.GHz]

for channel_id in dipole_channels:
    passband_low_per_channel[channel_id] = {}
    passband_high_per_channel[channel_id] = {}
    passband_low_per_channel[channel_id] = [0, 220 * units.MHz]
    passband_high_per_channel[channel_id] = [96 * units.MHz, 100 * units.GHz]

for channel_id in shallow_channels:
    passband_low_per_channel[channel_id] = {}
    passband_high_per_channel[channel_id] = {}
    passband_low_per_channel[channel_id] = [0 * units.MHz, 100 * units.GHz]
    passband_high_per_channel[channel_id] = [0, 1000 * units.MHz, 100 * units.GHz]

sampling_rate_pa = 0.5 * units.GHz
# x2 for upsampling
window_4ant = int(16 * units.ns * sampling_rate_pa * 2.0)
step_4ant = int(8 * units.ns * sampling_rate_pa * 2.0)
main_low_angle = np.deg2rad(-59.54968597864437)
main_high_angle = np.deg2rad(59.54968597864437)
phasing_angles_4ant = np.arcsin(np.linspace(np.sin(main_low_angle), np.sin(main_high_angle), 11))

if __name__ == "__main__":
    results_folder = '/lustre/fs22/group/radio/lpyras/muon_sim/M2_results'
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    parser = argparse.ArgumentParser(description='Run NuRadioMC simulation')
    parser.add_argument('--inputfilename', type=str, default='/lustre/fs22/group/radio/lpyras/muon_sim/M1_simulation_input/input_1.0e+18_1.0e+19_cos_theta_0.1_0.0_n_100_job_0.hdf5',
                        help='path to NuRadioMC input event list')
    parser.add_argument('--detectordescription', type=str, default='/afs/ifh.de/group/radio/scratch/lpyras/muon_sim'
                                                                   '/detector_4PA_single_station.json',
                        help='path to file containing the detector description')
    parser.add_argument('--config', type=str, default='/afs/ifh.de/group/radio/scratch/lpyras/muon_sim/config.yaml',
                        help='NuRadioMC yaml config file')
    parser.add_argument('--outputfilename', type=str, default=os.path.join(results_folder, 'MC_input_1.0e+18_1.0e+19_cos_theta_0.1_0.0_n_100_job_2.hdf5'),
                        help='hdf5 output filename')
    parser.add_argument('--outputfilename_NuRadioReco', type=str, nargs='?', default=None,
                        help='outputfilename of NuRadioReco detector sim file')
    args = parser.parse_args()

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
            channelBandPassFilter.run(evt, station_copy, det, passband=passband_low_per_channel, filter_type="cheby1", order=9, rp=.1)
            channelBandPassFilter.run(evt, station_copy, det, passband=passband_high_per_channel, filter_type="cheby1", order=4, rp=.1)

            Vrms_deep = np.sqrt(np.mean((station_copy.get_channel(deep_channels[0]).get_trace()) ** 2))
            Vrms_dipole = np.sqrt(np.mean((station_copy.get_channel(dipole_channels[0]).get_trace()) ** 2))
            Vrms_shallow = np.sqrt(np.mean((station_copy.get_channel(shallow_channels[0]).get_trace()) ** 2))

            print('Vrms_deep', Vrms_deep)
            print('Vrms_dipole', Vrms_dipole)
            print('Vrms_shallow', Vrms_shallow)

            # LPDAs
            highLowThreshold.run(evt, station_copy, det,
                                 threshold_high=2.5 * Vrms_shallow,
                                 threshold_low=-2.5 * Vrms_shallow,
                                 coinc_window=50 * units.ns,
                                 triggered_channels=shallow_channels,
                                 number_concidences=2,  # 2/3 majority logic
                                 trigger_name='hilo_2of3_2.5sigma')

            highLowThreshold.run(evt, station_copy, det,
                                 threshold_high=3 * Vrms_shallow,
                                 threshold_low=-3 * Vrms_shallow,
                                 coinc_window=50 * units.ns,
                                 triggered_channels=shallow_channels,
                                 number_concidences=2,  # 2/3 majority logic
                                 trigger_name='hilo_2of3_3sigma')

            highLowThreshold.run(evt, station_copy, det,
                                 threshold_high=3.5 * Vrms_shallow,
                                 threshold_low=-3.5 * Vrms_shallow,
                                 coinc_window=50 * units.ns,
                                 triggered_channels=shallow_channels,
                                 number_concidences=2,  # 2/3 majority logic
                                 trigger_name='hilo_2of3_3.5sigma')
            # noiseless dipole
            simpleThreshold.run(evt, station_copy, det,
                                threshold=1.5 * Vrms_dipole,
                                triggered_channels=dipole_channels,
                                trigger_name=f'dipole_1.5sigma')

            simpleThreshold.run(evt, station_copy, det,
                                threshold=2 * Vrms_dipole,
                                triggered_channels=dipole_channels,
                                trigger_name=f'dipole_2sigma')

            simpleThreshold.run(evt, station_copy, det,
                                threshold=2.5 * Vrms_dipole,
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

            # if(station_copy.has_triggered()):
            #     print("TRIGGERED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", flush=True)
            # 5) set trigger attributes of original station
            for trigger in station_copy.get_triggers().values():
                # print(trigger)
                station.set_trigger(trigger)

            # this module cuts the trace to the record length of the detector
            triggerTimeAdjuster.run(evt, station, det)


    sim = mySimulation(inputfilename=args.inputfilename,
                       outputfilename=args.outputfilename,
                       detectorfile=args.detectordescription,
                       outputfilenameNuRadioReco=args.outputfilename_NuRadioReco,
                       config_file=args.config)
    sim.run()
