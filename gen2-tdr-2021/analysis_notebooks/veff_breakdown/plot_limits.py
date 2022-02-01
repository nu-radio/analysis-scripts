import numpy as np
import matplotlib.pyplot as plt
from NuRadioMC.examples.Sensitivities import E2_fluxes3 as limits
from NuRadioMC.utilities import fluxes
from NuRadioReco.utilities import units
from NuRadioMC.utilities import cross_sections
import helper
import plot_helper as plot_helper

"""
This script takes the results of plot_veffs.py
and makes a limit plot. It also includes the optical EHE values.
Also, it will make an effective area plot.

"""


livetime = 10 * units.year
uptime = 0.9

detectors = [
    "review_array",
    "baseline_array",
    "hex_hybrid_only_array",
    "hex_shallow_array",
    "hex_shallowheavy_array",
]

config = "config_ARZ2020_noise"
detsim = "D01detector_sim"


fig, ax = limits.get_E2_limit_figure(show_ice_cube_EHE_limit=True,
										show_ice_cube_HESE_data=False,
										show_ice_cube_HESE_fit=False,
										show_ice_cube_mu=True,
										show_anita_I_IV_limit=True,
										show_auger_limit=False,
										diffuse=True, show_grand_10k=False, 
                                        show_grand_200k=False, show_Heinze=False, show_TA=False,
										show_ara=False, show_arianna=False, shower_Auger=False)

# for veff
fig2, ax_veff = plt.subplots(1, 2, figsize=(12,6))

# for aeff
fig3, ax_aeff = plt.subplots(1, 1)

linestyle=['dotted', 'dashed', 'dashdot', (0, (5, 10))]

########################
########################
# plot the review array
########################
########################

det = 'review_array'
deep_trigger = 'review'
shallow_trigger = 'review'

data_review = np.genfromtxt(f'results/veff_{det}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}_mode_review.csv', 
    delimiter=',', skip_header=1, 
    names=['logE', 'veff', 'dveff', 'sveff', 'coincveff', 'aeff', 'daeff', 'saeff', 'coincaeff', 'fdeep', 'fshallow', 'fcoinc'])
energies = np.power(10.,data_review['logE'])
limit_review = fluxes.get_limit_e2_flux(energy=energies*units.eV, 
                                        veff_sr=data_review['veff']*units.km**3*units.sr, 
                                        livetime=livetime*uptime, upperLimOnEvents=1)
ax.plot(energies/limits.plotUnitsEnergy, limit_review/limits.plotUnitsFlux, 
    'C4-', linewidth=3,
    label='{}, {} yrs, {}% uptime'.format(det,livetime/units.year, uptime*100))
ax_veff[0].plot(energies, data_review['veff'],
    'C4-',
    label='{}'.format(det)
    )

review_aeff = []
for e, v in zip(energies*units.eV, data_review['veff']*units.km**3*units.sr):
    lint =  cross_sections.get_interaction_length(float(e))
    aeff_temp = v/lint
    review_aeff.append(aeff_temp)
ax_aeff.plot(energies, review_aeff, 
    linewidth=2,
    label='Feb 2021 Radio Review Array, Trigger Level')

########################
########################
# now plot the new arrays
########################
########################

deep_trigger = 'PA_4channel_100Hz'
shallow_trigger = 'LPDA_2of4_100Hz'

for i, det in enumerate(detectors[1:]):
    data = np.genfromtxt(f'results/veff_{det}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}_mode_deepshallow.csv', 
        delimiter=',', skip_header=1, 
        names=['logE', 'veff', 'dveff', 'sveff', 'coincveff', 'aeff', 'daeff', 'saeff', 'coincaeff', 'fdeep', 'fshallow', 'fcoinc'])
    energies = np.power(10.,data['logE'])

    limit = fluxes.get_limit_e2_flux(energy=energies*units.eV, 
                                        veff_sr=data['veff']*units.km**3*units.sr, 
                                        livetime=livetime*uptime, upperLimOnEvents=1)
    l = ax.plot(energies/limits.plotUnitsEnergy, 
        limit/limits.plotUnitsFlux, 
        linewidth=3, 
        linestyle=linestyle[i],
        label='{}, {} yrs, {}% uptime'.format(det,livetime/units.year, uptime*100)
        )

    ax_veff[0].plot(energies, data['veff'],
        linestyle=linestyle[i], linewidth=3,
        label='{}'.format(det)
        )

    ax_veff[1].plot(energies, data['veff']/data_review['veff'],
        linestyle=linestyle[i],linewidth=3,
        label='{}'.format(det)
    )

########################
########################
# sketch in the optical array
########################
########################

optical_energies, optical_aeff = plot_helper.get_gen2opticalehe()
optical_energies = optical_energies*units.GeV
optical_aeff = optical_aeff*units.m**2 * units.sr
optical_veff = []
for e, aeff in zip(optical_energies[1:], optical_aeff):
    lint = cross_sections.get_interaction_length(float(e))

    the_veff = aeff * lint
    optical_veff.append(the_veff/units.km**3)

ax_aeff.plot(*plot_helper.stepped_path(optical_energies, optical_aeff),
    linewidth=2,
    label='Gen2-Optical EHE, Analysis Level'
    )

ax_veff[0].plot(*plot_helper.stepped_path(optical_energies, optical_veff))


limit_optical = fluxes.get_limit_e2_flux(energy=optical_energies[1:]/units.eV, 
                                        veff_sr=np.asarray(optical_veff)*units.km**3*units.sr,
                                        livetime=livetime*uptime, upperLimOnEvents=1)
ax.plot(optical_energies[1:]/limits.plotUnitsEnergy, limit_optical/limits.plotUnitsFlux, 
    'C5-', linewidth=3,
    label='{}, {} yrs, {}% uptime'.format('optical',livetime/units.year, uptime*100))



########################
########################
# science goals
########################
########################

ax.plot([np.power(10., 17.5), np.power(10., 19)], [3E-10, 3E-10], 'r-', linewidth=3)
ax.plot([30E15], [1.9E-9], 'rx', markersize=10, markeredgewidth=3)


########################
########################
# beautify and save
########################
########################

ax.legend(loc='upper right')
ax.set_ylim([1E-10, 1E-7])
ax.set_xlim([1E15, 1E21])
ax.set_title('Trigger Level')
fig.savefig('plots/limits.png')

ax_veff[0].set_xscale('log')
ax_veff[0].set_yscale('log')
ax_veff[0].set_title("Trigger Level Effective Volume")
ax_veff[0].set_xlabel('neutrino energy eV')
ax_veff[0].set_ylabel(r'Effective Volume [km$^3$ str w.e.]')
ax_veff[0].set_ylim([1E-1, 1E4])
ax_veff[0].set_xlim([1E14, 1E20])

ax_veff[1].set_xscale('log')
ax_veff[1].set_title("New Array Veff / Feb 2021 Review Array Veff")
ax_veff[1].set_xlabel('neutrino energy eV')
ax_veff[1].set_ylabel(r'Ratio')
ax_veff[1].legend(loc='upper right')

fig2.tight_layout()
fig2.savefig('plots/all_veffs.png')


# the effective areas plots
ax_aeff.legend()
ax_aeff.set_xlim([1E15, 1E20])
ax_aeff.set_ylim([1E1, 1E8])
ax_aeff.set_xscale('log')
ax_aeff.set_yscale('log')
ax_aeff.set_title("All-Sky, Flavor-Averaged Effecive Area")
ax_aeff.set_xlabel('neutrino energy eV')
ax_aeff.set_ylabel(r'Effective Area [m$^2$ str]')
plt.tight_layout()
fig3.savefig('plots/aeff_comparison.png')
