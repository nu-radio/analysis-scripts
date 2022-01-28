import numpy as np
import matplotlib.pyplot as plt
from NuRadioReco.utilities import units
import plot_helper as helper


# for aeff
fig3, ax_aeff = plt.subplots(1, 1)

########################
########################
# plot the review array
########################
########################

radio_energies, radio_veff, radio_aeff = helper.get_radio_review_array()
radio_energies = radio_energies*units.eV
radio_aeff = radio_aeff*units.km**2*units.sr

ax_aeff.plot(radio_energies/units.eV, radio_aeff/units.m**2/units.sr, 
    linewidth=2,
    label='Feb 2021 Radio Review Array, Trigger Level')

########################
########################
# sketch in the optical array
########################
########################

optical_energies, optical_aeff = helper.get_gen2opticalehe()
optical_energies = optical_energies*units.GeV
optical_aeff = optical_aeff*units.m**2 * units.sr

ax_aeff.plot(*helper.stepped_path(optical_energies/units.eV, optical_aeff/units.m**2/units.sr),
    linewidth=2,
    label='Gen2-Optical EHE, Analysis Level'
    )


########################
########################
# beautify and save
########################
########################

ax_aeff.legend()
ax_aeff.set_xlim([1E15, 1E20])
ax_aeff.set_ylim([1E1, 1E8])
ax_aeff.set_xscale('log')
ax_aeff.set_yscale('log')
ax_aeff.set_title("All-Sky, Flavor-Averaged Effecive Area")
ax_aeff.set_xlabel('Neutrino Energy (eV)')
ax_aeff.set_ylabel(r'Effective Area [m$^2$ str]')
plt.tight_layout()
fig3.savefig('aeff_comparison.png')
