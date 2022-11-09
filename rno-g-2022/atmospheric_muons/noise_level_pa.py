import argparse
import time
import numpy as np
from astropy.time import Time
from multiprocessing import Pool as ThreadPool
import NuRadioReco.modules.channelGenericNoiseAdder
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.phasedarray.triggerSimulator
import NuRadioReco.modules.trigger.simpleThreshold
from NuRadioReco.framework.event import Event
from NuRadioReco.framework.station import Station
from NuRadioReco.framework.channel import Channel
from NuRadioReco.utilities import units
from NuRadioReco.utilities import fft
from NuRadioReco.detector import detector
import NuRadioReco.modules.RNO_G.hardwareResponseIncorporator


parser = argparse.ArgumentParser(description='calculates noise trigger rate for phased array')
parser.add_argument('--ntries', type=int, help='number noise traces to which a trigger is applied for each threshold',
                    default=1000)
parser.add_argument('--ncpus', type=int, help='number of parallel jobs that can be run',
                    default=1)
parser.add_argument('--detectordescription', type=str, default= 'detector_PA.json',
                    help='path to file containing the detector description')
parser.add_argument('--nchannels', type=int,
                    help='number of channels to phase', default=4)

args = parser.parse_args()


thresholds = np.arange(35, 45, 1)


det_file = args.detectordescription
n_channels = args.nchannels
ntrials = args.ntries
ncpus = args.ncpus

det = detector.Detector(json_filename=det_file)
det.update(Time.now())


main_low_angle = np.deg2rad(-59.55)
main_high_angle = np.deg2rad(59.55)

station_id = 11
channel_ids = det.get_channel_ids(station_id)

n_samples = det.get_number_of_samples(station_id, channel_ids[0])  # we assume that all channels have the same parameters
sampling_rate = det.get_sampling_frequency(station_id, channel_ids[0])
dt = 1 / sampling_rate

print(f"sampling rate = {sampling_rate/units.MHz}MHz, {n_samples} samples")

if(n_channels == 4):
    upsampling_factor = 2
    default_angles = np.arcsin(np.linspace(np.sin(main_low_angle), np.sin(main_high_angle), 11))
else:
    print("wrong n_channels!")
    exit()

window_length = int(16 * units.ns * sampling_rate * upsampling_factor)
step_length = int(8 * units.ns * sampling_rate * upsampling_factor)

pattern = f"pa_trigger_rate_{n_channels:d}channels_{upsampling_factor}xupsampling"

hardwareResponseIncorporator = NuRadioReco.modules.RNO_G.hardwareResponseIncorporator.hardwareResponseIncorporator()
triggerSimulator = NuRadioReco.modules.phasedarray.triggerSimulator.triggerSimulator()
thresholdSimulator = NuRadioReco.modules.trigger.simpleThreshold.triggerSimulator()
channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()


def loop(zipped):

    threshold = float(zipped[0])
    seed = int(zipped[1])

    evt = Event(0, 0)
    station = Station(station_id)
    channel_ids = det.get_channel_ids(station_id)
    for channel_id in channel_ids:  # take some channel id that match your detector
        channel = Channel(channel_id)
        default_trace = np.zeros(8192)
        channel.set_trace(trace=default_trace, sampling_rate=sampling_rate)
        station.add_channel(channel)

    channelGenericNoiseAdder = NuRadioReco.modules.channelGenericNoiseAdder.channelGenericNoiseAdder()
    channelGenericNoiseAdder.begin(seed=seed)

    # Noise rms is amplified to greater than Vrms so that, after filtering, its the Vrms we expect
    channelGenericNoiseAdder.run(evt, station, det,
                             amplitude=1.2211064327076491e-05,
                             min_freq=80*units.MHz,
                             max_freq=800*units.MHz, type='rayleigh', bandwidth=None)
    
    hardwareResponseIncorporator.run(evt, station, det, sim_to_data=True)
    
    #filter
    channelBandPassFilter.run(evt, station, det, passband=[0, 220 * units.MHz], filter_type="cheby1", order=9, rp=.1)
    channelBandPassFilter.run(evt, station, det, passband=[96 * units.MHz, 100 * units.GHz], filter_type="cheby1", order=4, rp=.1)
    
    trace = (station.get_channel(0).get_trace())
    Vrms = np.sqrt(np.mean(trace**2))
    threshold_ = threshold * np.power(Vrms, 2.0)

    triggered = triggerSimulator.run(evt, station, det,
                                     Vrms,
                                     threshold_,
                                     triggered_channels=channel_ids,
                                     phasing_angles=default_angles,
                                     ref_index=1.76,
                                     trigger_name='primary_phasing',
                                     trigger_adc=False,
                                     adc_output='voltage',
                                     trigger_filter=None,
                                     upsampling_factor=upsampling_factor,
                                     window=window_length,
                                     step=step_length)

    return triggered


pool = ThreadPool(ncpus)

for threshold in thresholds:
    print(thresholds)
    print('tested threshold', threshold * np.power(0.005526947222857975, 2.0))
    n_triggers = 0
    i = 0
    t00 = time.time()
    t0 = time.time()

    while i < ntrials:
        n_pool = int(float(ntrials) / 10.0)

        print("Events:", i, " Delta t=", time.time() - t00, " N_triggers =", n_triggers)
        i += n_pool
        t00 = time.time()

        results = pool.map(loop, zip(threshold * np.ones(n_pool), np.random.get_state()[1][0] + i + np.arange(n_pool)))
        n_triggers += np.sum(results)
        rate = 1. * n_triggers / (i * n_samples * dt)

    rate = 1. * n_triggers / (i * n_samples * dt)

    with open(f"{pattern}.txt", "a") as fout:
        fout.write(f"{threshold}\t{n_triggers}\t{i*n_samples*dt}\t{rate}\n")
        fout.close()
    print(f"threshold = {threshold:.3f}: n_triggers = {n_triggers} -> rate = {rate/units.Hz:.0f} Hz, {(time.time() -  t0)/i*1000:.1f}ms per event -> {(time.time() -  t0)/n_triggers*100/60:.1f}min for 100 triggered events")
