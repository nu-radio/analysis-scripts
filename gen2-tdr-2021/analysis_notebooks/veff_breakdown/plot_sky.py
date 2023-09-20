import os
import numpy as np
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

top_dir = '/data/sim/Gen2/radio/2020/gen2-tdr-2021/simulation_output/secondaries_1700km2/step4/fine_zenithbinning'
arrays = [
    'baseline_array',
    'hex_hybrid_only_array',
    'hex_shallow_array',
    'hex_shallowheavy_array'
]

thetamins = None

results = {}
for a in arrays:
    results[a] = np.zeros(80)
energy_choice = 1E18

fnames = "{}/*.json".format(top_dir)
files = glob(fnames)
dataframes = []
for f in files:
    df = pd.read_json(f)
    a, detector, c, d, flavor = f.split("__")
    flavor = flavor.split('.')[0]
    print("Ingesting Det {}, Flavor {}".format(detector, flavor))
    df.sort_values(["thetamin"], inplace=True) # sort by thetamin
    df.reset_index(inplace=True) # reset the row index
    energy_mask = df['energy']==energy_choice # grab only the 1 EeV slice
    selection = df.loc[energy_mask]

    thetamin = selection.thetamin
    veff = selection.veff

    if thetamins is None:
        thetamins = thetamin

    tempy = np.zeros(80)

    for i, v in enumerate(veff):
        tempy[i]=v['combined_4channelPA'][0]*1E-9 # convert to km^3
    
    results[detector] += tempy/3.

fig, ax = plt.subplots(1, 1)
linestyles=['-', '-.', '--', '-.']
linecolors=['C0', 'C1', 'C2', 'C3']

# arrays = [
#     'baseline_array'
# ]


for i, a in enumerate(arrays):

    the_y_values = results[a]
    max_value = np.max(the_y_values)/2.
    the_cut = the_y_values > max_value

    coszen_values = np.cos(thetamins)
    sindec_values = -coszen_values
    bin_width = abs((coszen_values.values[1] - coszen_values.values[0])/2.)
    min_coszen_val = np.min(coszen_values[the_cut])-bin_width
    max_coszen_val = np.max(coszen_values[the_cut])+bin_width
    max_sindec_val = -min_coszen_val
    min_sindec_val = -max_coszen_val
    # print("Array {}, Min CosZen {}, Max CosZen {}, MinSinDec {}, MaxSinDec {}".format(a,min_coszen_val, max_coszen_val, min_sindec_val, max_sindec_val))
    # ax.axvline(min_sindec_val, color=linecolors[i], label='Min Sin Dec')
    # ax.axvline(max_sindec_val, color=linecolors[i], label='Max Sin Dec')
    sky_area = (180./np.pi) * 360. * (max_sindec_val - min_sindec_val)
    sky_fraction = sky_area / 41252.96125 # fraction of whole sky
    print("Array {}, Sky Area {}, Sky Fraction {}".format(a, sky_area, sky_fraction))

    ax.plot(
        sindec_values,
        results[a],
        linestyle=linestyles[i],
        color=linecolors[i],
        label='{}'.format(a)
    )


# plot the galactic center
dec_gc = 29.
print("Sin Dec of GC {}".format(np.sin(dec_gc)))
ax.axvline(np.sin(dec_gc), color='grey', linestyle='--')

ax.legend()
ax.set_xlabel(r'sin($\delta$)', fontsize=16)
ax.set_ylabel(r'V$_{eff}$ [km$^3$]', fontsize=16)
ax.set_yscale('log')
ax.set_ylim([5E1, 1.1E3])
# ax.set_xlim([-0.25, 1])
ax.set_xlim([-1, 0.25])
ax.set_title("Energy: {}".format(energy_choice))

plt.tight_layout()
fig.savefig('zen_plots.png')
del fig, ax


    