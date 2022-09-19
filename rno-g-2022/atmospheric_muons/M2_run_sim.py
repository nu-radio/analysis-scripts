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

# deep config
deep_channels = [0, 1, 2, 3]
threshold_pa = 38.62
sampling_rate_pa = 0.5 * units.GHz
# run the 4 phased trigger
# x2 for upsampling
window_4ant = int(16 * units.ns * sampling_rate_pa * 2.0)
step_4ant = int(8 * units.ns * sampling_rate_pa * 2.0)

main_low_angle = np.deg2rad(-59.54968597864437)
main_high_angle = np.deg2rad(59.54968597864437)
phasing_angles_4ant = np.arcsin(np.linspace(np.sin(main_low_angle), np.sin(main_high_angle), 11))

shallow_channels = [12, 14, 15, 17, 18, 20]
thresholds_LPDA = {
    '2/4_100Hz': 3.9498194908011524,
    '2/4_10mHz': 4.919151494949084}


if __name__ == "__main__":
    results_folder = '/lustre/fs22/group/radio/lpyras/muon_sim/M2_results'
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    parser = argparse.ArgumentParser(description='Run NuRadioMC simulation')
    parser.add_argument('--inputfilename', type=str, default='/lustre/fs22/group/radio/lpyras/muon_sim'
                                                      '/M1_simulation_input/input_1.0e+19_1.0e+20_cos_theta_1.0_0.0_job_0.hdf5',
                        help='path to NuRadioMC input event list')
    parser.add_argument('--detectordescription', type=str, default='/afs/ifh.de/group/radio/scratch/lpyras/muon_sim'
                                                                   '/detector.json',
                        help='path to file containing the detector description')
    parser.add_argument('--config', type=str, default='/afs/ifh.de/group/radio/scratch/lpyras/muon_sim/config.yaml',
                        help='NuRadioMC yaml config file')
    parser.add_argument('--outputfilename', type=str, default=os.path.join(results_folder, 'MC_input_1.0e+19_1.0e+20_cos_theta_1.0_0.0_job_0.hdf5'),
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
            #highLowThreshold.run(evt, station, det,
            #                     threshold_high=2 * self._Vrms,
            #                     threshold_low=-2 * self._Vrms,
            #                     coinc_window=40 * units.ns,
            #                     triggered_channels=shallow_channels,
            #                     number_concidences=2,  # 2/4 majority logic
            #                     trigger_name='hilo_2of4_2_sigma')

            simpleThreshold.run(evt, station, det,
                                threshold=.5 * self._Vrms,
                                triggered_channels=shallow_channels,
                                trigger_name='simple_05_sigma')

            phasedArrayTrigger.run(evt, station, det,
                                   Vrms=self._Vrms,
                                   threshold=threshold_pa * np.power(self._Vrms, 2.0),
                                   triggered_channels=deep_channels,
                                   phasing_angles=phasing_angles_4ant,
                                   ref_index=1.75,
                                   trigger_name=f'PA_4channel_100Hz',  # the name of the trigger
                                   trigger_adc=True,  # Don't have a seperate ADC for the trigger
                                   adc_output=f'voltage',  # output in volts
                                   trigger_filter=None,
                                   upsampling_factor=2,
                                   window=window_4ant,
                                   step=step_4ant)


    sim = mySimulation(inputfilename=args.inputfilename,
                       outputfilename=args.outputfilename,
                       detectorfile=args.detectordescription,
                       outputfilenameNuRadioReco=args.outputfilename_NuRadioReco,
                       config_file=args.config,
                       default_detector_station=11,
                       default_detector_channel=1)
    sim.run()
