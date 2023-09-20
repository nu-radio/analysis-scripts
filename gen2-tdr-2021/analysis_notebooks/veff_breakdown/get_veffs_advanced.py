import fractions
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import h5py
import glob
import json
import os
import sys
import pickle
from NuRadioReco.utilities import units
from NuRadioMC.utilities import cross_sections
import helper
import argparse
from radiotools import plthelpers as php


"""
This script caluculates takes the output of the calc_fractions_advanced.py,
and makes (1) a csv file of the resultant veffs and fractions,
as well as a plot of the Veff and fractions.

There are currently two supported modes:
stations_and_triggers and triggers_coincidences
Read the helper for a more thorough description of what each one does 

"""

parser = argparse.ArgumentParser()
parser.add_argument('--det', required=True, help='which detector index')
parser.add_argument('--mode', required=True, help='which mode: deepshallow or hybridshallow')
args = parser.parse_args()

the_det = int(args.det)
mode = args.mode
helper.check_requested_mode(mode)
the_list, labels_dict = helper.get_list_and_labels(mode)

detectors = helper.get_detectors()
flavors = helper.get_flavors()

top = "/data/sim/Gen2/radio/2020/gen2-tdr-2021/simulation_output/secondaries_1700km2/step3"
detector = detectors[the_det]
config = "config_ARZ2020_noise"
detsim = "D01detector_sim"
deep_trigger = 'PA_4channel_100Hz'
shallow_trigger = 'LPDA_2of4_100Hz'


# calculate the result per flavor
result = {}
for flavor in flavors:
    result[flavor] = {}
    result[flavor]['E'] = []
    for l in the_list:
        result[flavor][l+'_veff'] = []
        result[flavor][l+'_aeff'] = []
    
    pkl_file_name = os.path.join('data/overlap_' + detector + "_deeptrig_" + deep_trigger + "_shallowtrig_" + shallow_trigger + "_mode_"+ mode+'_' + flavor + ".pkl")
    with open(pkl_file_name, "rb") as fin:
        coincidences = pickle.load(fin)
    
    for lgE, value in coincidences.items():
        result[flavor]['E'].append(lgE)
        lint = cross_sections.get_interaction_length(np.power(10., float(lgE)))

        for l in the_list:
            v_temp = value[l] # the veff water equivalent
            a_temp = v_temp/lint # get the aeff water equivalent
            v_temp = v_temp * 4 * np.pi / units.km**3 # convert to km^3 steradians
            a_temp = a_temp * 4 * np.pi / units.km**2 # convert to km^2 steradians

            result[flavor][l+'_veff'].append(v_temp)
            result[flavor][l+'_aeff'].append(a_temp)

# average over flavors
energies = result['e']['E']
result_average = {}
for l in the_list:
    result_average[l+'_veff'] = np.zeros(9)
    result_average[l+'_aeff'] = np.zeros(9)
for l in the_list:
    for iflavor, flavor in enumerate(flavors):
        for iE, energy in enumerate(energies):
            result_average[l+'_veff'][iE] += result[flavor][l+'_veff'][iE]/3
            result_average[l+'_aeff'][iE] += result[flavor][l+'_aeff'][iE]/3

# calculate averages
result_fractions = {}
for l in the_list:
    result_fractions[l] = result_average[l+'_veff']/result_average['total_veff']

######################
# Plot Effective Area vs Energy
# And fractiosn vs Energy
######################

fig, axs = plt.subplots(1, 2, figsize=(12,6))
xx = result[flavor]['E']
fig.suptitle(f"{detector}")

for i, l in enumerate(the_list):
    axs[0].plot(xx, result_average[l+'_aeff'], php.get_color_linestyle(i), label=l)
    axs[1].plot(xx, result_fractions[l])

axs[0].set_yscale('log')
axs[0].set_xlabel("log10(energy [eV])")
axs[0].set_ylabel(r"[km$^2$ * str]")
axs[0].set_title("Water Equivalent Effective Area ")
axs[0].set_ylim([1E-5,1E2])
axs[0].legend(loc='lower right')

axs[1].set_ylim([0, 1.1])
axs[1].set_xlabel("log10(energy [eV])")
axs[1].set_ylabel("Fraction")
axs[1].set_title("Fraction of the All-Sky Effective Area")
fig.tight_layout(pad=2)
outfilename = f'plots/veff_fractions_{detector}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}_mode_{mode}.png'
fig.savefig(outfilename)

del fig, axs



######################
# plot fraction vs Energy (bigger standalone plot)
######################

import matplotlib as mpl
plt.rcParams.update({
    'xtick.labelsize': 14, 
    'ytick.labelsize': 14,
    'xtick.major.size': 5, 
    'ytick.major.size': 5,
    'axes.titlesize': 14,
    'axes.labelsize': 14
})
mpl.rc('font', size=12)
mpl.rc('axes', titlesize=14)

fig, axs = plt.subplots(1, 1, figsize=(12,6))
for i, l in enumerate(the_list):
    axs.plot(xx, result_fractions[l], php.get_color_linestyle(i), label=labels_dict[l])

axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs.set_ylim([0, 1.1])
axs.set_xlabel("log10(energy [eV])")
axs.set_ylabel("Fraction")
axs.set_title("Fraction of the All-Sky Effective Area")
fig.suptitle(f"{detector}")
fig.tight_layout()
outfilename = f'plots/fractions_{detector}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}_mode_{mode}.png'
fig.savefig(outfilename)



######################
# Save to CSV File
######################

# build the header file
output_csv_3 = "energy, veff [km^3]"
for l in the_list[1:]:
    label = labels_dict[l]
    output_csv_3+=", {}".format(label)
output_csv_3+="\n"

for iE, energy in enumerate(energies):
    tmp = ''
    tmp+='{:.1f}'.format(float(energy))
    tmp+=', {:e}'.format(result_average['total_veff'][iE])
    for l in the_list[1:]:
        tmp+=', {:>.4f}'.format(result_fractions[l][iE])
    tmp+='\n'
    output_csv_3+=tmp
print(output_csv_3)

# # now, build the formatting string itself
# formatting_string = '{:.1f}, {:e}' # always start with the energy and overal veff in scientific notation
# for l in the_list[1:]:
#     addition = ', {:>.4f}'
#     formatting_string+=addition


with open(f'results/veff_fractions_{detector}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}_{mode}.csv', 'w') as fout:
    fout.write(output_csv_3)


