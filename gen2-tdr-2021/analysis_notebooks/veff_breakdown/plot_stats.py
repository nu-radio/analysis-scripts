import numpy as np
import matplotlib.pyplot as plt

"""
This script takes the results of plot_veffs.py
and makes histograms of the number of events thrown, 
and number of events triggered.

Run like: python plot_stats.py

"""

detectors = [
    "baseline_array",
    "hex_hybrid_only_array",
    "hex_shallow_array",
    "hex_shallowheavy_array",
]
deep_trigger = 'PA_4channel_100Hz'
shallow_trigger = 'LPDA_2of4_100Hz'

for det in detectors:
    data = np.genfromtxt(f'results/stats_{det}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}.csv', 
        delimiter=',', skip_header=1, 
        names=['logE', 'n_thrown', 'trig_weight'])
    energies = np.power(10.,data['logE'])

    strings = []
    for e in data['logE']:
        strings.append("{}".format(e))

    fig, ax = plt.subplots(1, 1)
    ax.bar(strings, data['n_thrown'], label='number thrown: {:.1e}'.format(np.sum(data['n_thrown'])))
    ax.bar(strings, data['trig_weight'], label='trigering weight: {:.1e}'.format(np.sum(data['trig_weight'])))


    ax.legend()
    ax.set_yscale('log')
    ax.set_title("Simulation Statistics, {}".format(det))
    ax.set_xlabel(r'log$_{10}(E_{\nu})$')
    ax.set_ylabel(r'Number of Events')
    plt.tight_layout()
    fig.savefig('plots/sim_stats_{}.png'.format(det))
    del fig, ax



