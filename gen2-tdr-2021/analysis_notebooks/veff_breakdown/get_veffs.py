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
import plot_helper as plot_helper
import argparse

"""
This script caluculates takes the output of the calc_fractions.py,
and makes (1) a csv file of the resultant veffs and fractions,
as well as a plot of the Veff and fractions.

"""


parser = argparse.ArgumentParser()
parser.add_argument('--det', required=True, help='which detector index')
parser.add_argument('--mode', required=True, help='which mode: deepshallow or hybridshallow')
args = parser.parse_args()

the_det = int(args.det)
the_mode = args.mode

detectors = [
    "baseline_array",
    "hex_hybrid_only_array",
    "hex_shallow_array",
    "hex_shallowheavy_array",
    "review_array"
]

top = "/data/sim/Gen2/radio/2020/gen2-tdr-2021/simulation_output/secondaries_1700km2/step3"
detector = detectors[the_det]
config = "config_ARZ2020_noise"
detsim = "D01detector_sim"
deep_trigger = 'PA_4channel_100Hz'
shallow_trigger = 'LPDA_2of4_100Hz'
mode = the_mode

do_review = False
if the_det==4:
    do_review = True

flavors = [
    "e",
    "mu",
    "tau"
]

energies = np.arange(16, 20.5, 0.5)
average_total_veff = np.zeros(9)
average_deeporhybrid_veff = np.zeros(9)
average_shallow_veff = np.zeros(9)
average_dual_veff = np.zeros(9)
average_total_aeff = np.zeros(9)
average_deeporhybrid_aeff = np.zeros(9)
average_shallow_aeff = np.zeros(9)
average_dual_aeff = np.zeros(9)
total_n_thrown = np.zeros(9)
total_trig_weight = np.zeros(9)

if not do_review:
    result = {}
    for flavor in flavors:
        result[flavor] = {
                        'E' : [],
                        'total_veff': [],
                        'deeporhybrid_veff': [],
                        'shallow_veff': [],
                        'dual_veff': [],
                        'total_aeff': [],
                        'deeporhybrid_aeff': [],
                        'shallow_aeff': [],
                        'dual_aeff': [],
                        'n_thrown': [],
                        'trig_weight': []
                        }
        
        pkl_file_name = os.path.join('data/overlap_' + detector + "_deeptrig_" + deep_trigger + "_shallowtrig_" + shallow_trigger + "_mode_"+ mode+'_' + flavor + ".pkl")
        with open(pkl_file_name, "rb") as fin:
            coincidences = pickle.load(fin)
        
        for lgE, value in coincidences.items():
            result[flavor]['E'].append(lgE)
            veff_total = value['total_veff']
            veff_deeporhybrid = value['deeporhybrid_only_veff']
            veff_shallow = value['shallow_only_veff']
            veff_dual = value['dual_veff']
            
            lint = cross_sections.get_interaction_length(np.power(10., float(lgE)))
            aeff_total = veff_total/lint
            aeff_deeporhybrid = veff_deeporhybrid/lint
            aeff_shallow = veff_shallow/lint
            aeff_dual = veff_dual/lint

            # convert to the right units (I could be more compact than this, but want to be explicit)
            veff_total = veff_total * 4 * np.pi / units.km**3
            veff_deeporhybrid = veff_deeporhybrid * 4 * np.pi / units.km**3
            veff_shallow = veff_shallow * 4 * np.pi / units.km**3
            veff_dual = veff_dual * 4 * np.pi / units.km**3

            aeff_total = aeff_total * 4 * np.pi / units.km**2
            aeff_deeporhybrid = aeff_deeporhybrid * 4 * np.pi / units.km**2
            aeff_shallow = aeff_shallow * 4 * np.pi / units.km**2
            aeff_dual = aeff_dual * 4 * np.pi / units.km**2

            result[flavor]['total_veff'].append(veff_total)
            result[flavor]['deeporhybrid_veff'].append(veff_deeporhybrid)
            result[flavor]['shallow_veff'].append(veff_shallow)
            result[flavor]['dual_veff'].append(veff_dual)
            result[flavor]['total_aeff'].append(aeff_total)
            result[flavor]['deeporhybrid_aeff'].append(aeff_deeporhybrid)
            result[flavor]['shallow_aeff'].append(aeff_shallow)
            result[flavor]['dual_aeff'].append(aeff_dual)
            result[flavor]['n_thrown'].append(value['n_thrown'])
            result[flavor]['trig_weight'].append(value['trig_weight'])


    for iflavor, flavor in enumerate(flavors):
        for iE, energy in enumerate(energies):
            average_total_veff[iE]+=result[flavor]['total_veff'][iE]/3
            average_deeporhybrid_veff[iE]+=result[flavor]['deeporhybrid_veff'][iE]/3
            average_shallow_veff[iE]+=result[flavor]['shallow_veff'][iE]/3
            average_dual_veff[iE]+=result[flavor]['dual_veff'][iE]/3
            average_total_aeff[iE]+=result[flavor]['total_aeff'][iE]/3
            average_deeporhybrid_aeff[iE]+=result[flavor]['deeporhybrid_aeff'][iE]/3
            average_shallow_aeff[iE]+=result[flavor]['shallow_aeff'][iE]/3
            average_dual_aeff[iE]+=result[flavor]['dual_aeff'][iE]/3
            total_n_thrown[iE]+=result[flavor]['n_thrown'][iE]
            total_trig_weight[iE]+=result[flavor]['trig_weight'][iE]

scaling = 1
if do_review:
    average_total_veff, average_deeporhybrid_veff, average_shallow_veff, average_dual_veff, average_total_aeff, average_deeporhybrid_aeff, average_shallow_aeff, average_dual_aeff = plot_helper.get_review_array_padded()
    detector = 'review_array'
    deep_trigger='review'
    shallow_trigger='review'
    scaling = 4 * np.pi

fraction_deeporhybrid = average_deeporhybrid_veff/average_total_veff
fraction_shallow = average_shallow_veff/average_total_veff
fraction_dual = average_dual_veff/average_total_veff

if mode == 'deepshallow':
    sub1 = 'Deep'
    sub2 = 'Shallow'
elif mode == 'hybridshallow':
    sub1 = 'Hybrid Station'
    sub2 = 'Surface Stations'
elif do_review:
    sub1='Review'
    sub2='Review'

output_csv = f'log10(energy [eV]), veff total [km^3 sr water equiv], veff {sub1} only, veff shallow only, veff coincidence, '
output_csv += f'aeff total [km^2 sr water equiv], aeff {sub1} only, aeff shallow only, aeff coincidence, '
output_csv += f'frac {sub1} only, frac shallow only, frac coincidence'
output_csv += "\n"
for iE, energy in enumerate(energies):
    output_csv += '{:.1f}, {:e}, {:e}, {:e}, {:e}, {:e}, {:e}, {:e}, {:e}, {:>.2f}, {:>.2f}, {:>.2f} \n'.format(
        float(energy), 
        average_total_veff[iE]/scaling, average_deeporhybrid_veff[iE], average_shallow_veff[iE], average_dual_veff[iE],
        average_total_aeff[iE], average_deeporhybrid_aeff[iE], average_shallow_aeff[iE], average_dual_aeff[iE],
        fraction_deeporhybrid[iE], fraction_shallow[iE], fraction_dual[iE])
with open(f'results/veff_{detector}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}_mode_{mode}.csv', 'w') as fout:
    fout.write(output_csv)

output_csv_2 = 'energy, num thrown, total trig weight \n'
for iE, energy in enumerate(energies):
    output_csv_2 += '{:.1f}, {}, {} \n'.format(
        float(energy), total_n_thrown[iE], total_trig_weight[iE])
with open(f'results/stats_{detector}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}.csv', 'w') as fout2:
    fout2.write(output_csv_2)




fig, axs = plt.subplots(1, 2, figsize=(12,6))
colors = ['C0', 'C1', 'C2']
markers = ['o', 's', '^']
styles = ['C0o-', 'C1s--', 'C2^-.', 'C3v:', 'C4>-', 'C5<--']
xx = energies
fig.suptitle(f"{detector}")

axs[0].plot(xx, average_total_aeff, styles[3], label='Sum')
axs[0].plot(xx, average_deeporhybrid_aeff, styles[0], label=f'{sub1}-only component')
axs[0].plot(xx, average_shallow_aeff, styles[1], label=f'{sub2}-only component')
axs[0].plot(xx, average_dual_aeff, styles[2], label=f'{sub1} + {sub2} coincidence')
axs[0].plot(xx, average_deeporhybrid_aeff + average_dual_aeff, styles[4], label=f'{sub1}-inclusive component')
axs[0].plot(xx, average_shallow_aeff + average_dual_aeff, styles[5], label=f'{sub2}-inclusive component')
axs[0].set_yscale('log')
axs[0].set_xlabel("log10(energy [eV])")
axs[0].set_ylabel(r"[km$^2$ * str]")
axs[0].set_title("Water Equivalent Effective Area ")
axs[0].set_ylim([1E-5,1E2])
axs[0].legend(loc='lower right')


axs[1].plot(xx, fraction_deeporhybrid, styles[0], label=f'{sub1}-only component')
axs[1].plot(xx, fraction_shallow, styles[1], label=f'{sub2}-only component')
axs[1].plot(xx, fraction_dual, styles[2], label=f'{sub1} + {sub2} coincidence')
axs[1].plot(xx, fraction_deeporhybrid+fraction_shallow+fraction_dual, styles[3], label='Sum')
axs[1].plot(xx, fraction_deeporhybrid+fraction_dual, styles[4], label=f'{sub1}-inclusive component')
axs[1].plot(xx, fraction_shallow+fraction_dual, styles[5], label='{sub2}-inclusive component')
axs[1].set_ylim([0, 1.1])
axs[1].set_xlabel("log10(energy [eV])")
axs[1].set_ylabel("Fraction")
axs[1].set_title("Fraction of the All-Sky Effective Area")
fig.tight_layout(pad=2)
outfilename = f'plots/veff_fractions_{detector}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}_mode_{mode}.png'
fig.savefig(outfilename)

## for christian
from radiotools import plthelpers as php
fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 2]})

fixer = 4*np.pi
average_deeporhybrid_veff/=fixer
average_shallow_veff/=fixer
average_dual_veff/=fixer

xx_new = []
for x in xx:
    xx_new.append(float(x))
xx = xx_new

xx = np.power(10., np.asarray(xx))
ax.plot(xx, average_total_veff/scaling, php.get_color_linestyle(0), label='Sum')
ax.plot(xx, average_deeporhybrid_veff, php.get_color_linestyle(1), label=f'{sub1}-only')
ax.plot(xx, average_shallow_veff, php.get_color_linestyle(2), label=f'{sub2}-only')
ax.plot(xx, average_dual_veff, php.get_color_linestyle(3), label=f'{sub1} + {sub2} coincidence')
ax.plot(xx, average_deeporhybrid_veff + average_dual_veff, php.get_color_linestyle(4), label=f'{sub1}-inclusive')
ax.plot(xx, average_shallow_veff + average_dual_veff, php.get_color_linestyle(5), label=f'{sub2}-inclusive')

ax2.plot(xx, fraction_deeporhybrid, php.get_color_linestyle(1), label=f'{sub1}-only component')
ax2.plot(xx, fraction_shallow, php.get_color_linestyle(2), label=f'{sub2}-only component')
ax2.plot(xx, fraction_dual, php.get_color_linestyle(3), label=f'{sub1} + {sub2} coincidence')
ax2.plot(xx, fraction_deeporhybrid+fraction_dual, php.get_color_linestyle(4), label=f'{sub1}-inclusive component')
ax2.plot(xx, fraction_shallow+fraction_dual, php.get_color_linestyle(5), label=f'{sub2}-inclusive component')

ax.semilogx(True)
ax.set_title(detector)
ax.semilogy(True)
ax.set_ylim(1E-1, 1E3)
ax2.set_xlabel(f"energy [eV]")
ax2.set_ylabel(f"contrib. to\ntotal Veff")
ax2.set_ylim(0, 1)
# ax2.legend()
ax.set_ylabel(f"effective volume [km^3]")
ax.set_xlim(1e16, 1e20)
ax.legend(fontsize='small')
fig.tight_layout()
fig.savefig(f"plots/christian_veff_energy_{detector}.png")
# plt.show()
# a =1/0
plt.close("all")

if mode == 'deepshallow':
    output_csv_3 = "energy, veff [km^3], deep inclusive, shallow inclusive, deep exclusive, shallow exclusive, coincidences"
    output_csv_3 += "\n"
    for iE, energy in enumerate(energies):
        output_csv_3 += '{:.1f}, {:e}, {:>.4f}, {:>.4f}, {:>.4f}, {:>.4f}, {:>.4f} \n'.format(
            float(energy), 
            average_total_veff[iE]/scaling, 
            fraction_deeporhybrid[iE]+fraction_dual[iE],
            fraction_shallow[iE]+fraction_dual[iE],
            fraction_deeporhybrid[iE],
            fraction_shallow[iE],
            fraction_dual[iE]
            )
    with open(f'results/veff_fractions_{detector}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}.csv', 'w') as fout:
        fout.write(output_csv_3)
