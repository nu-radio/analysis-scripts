import numpy as np
import matplotlib.pyplot as plt
from NuRadioReco.utilities import units
import plot_helper as helper

"""
This script compares the Feb 2021 baseline effective area,
to that of the Gen2 Extremely High Energy (EHE)
effective area at analysis level.

So, be careful: this will compare radio trigger level
to optical analysis level.

"""


fig, ax = plt.subplots(1, 1)

########################
########################
# plot the review array
########################
########################

radio_energies, radio_veff, radio_aeff = helper.get_radio_review_array()
radio_energies = radio_energies*units.eV
radio_aeff = radio_aeff*units.km**2*units.sr

ax.plot(radio_energies/units.eV, radio_aeff/units.m**2/units.sr, 
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

ax.plot(*helper.stepped_path(optical_energies/units.eV, optical_aeff/units.m**2/units.sr),
    linewidth=2,
    label='Gen2-Optical EHE, Analysis Level'
    )


########################
########################
# beautify and save
########################
########################

ax.legend()
ax.set_xlim([1E15, 1E20])
ax.set_ylim([1E1, 1E8])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title("All-Sky, Flavor-Averaged Effecive Area")
ax.set_xlabel('Neutrino Energy (eV)')
ax.set_ylabel(r'Effective Area [m$^2$ str]')
plt.tight_layout()
fig.savefig('aeff_comparison.png')
