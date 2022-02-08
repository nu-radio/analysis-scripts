import sys
from NuRadioReco.utilities import units
from NuRadioReco.detector import detector

from NuRadioMC.utilities import fluxes
from NuRadioMC.utilities.Veff import  get_Veff_Aeff, get_Veff_Aeff_array, get_index, get_Veff_water_equivalent
from NuRadioMC.examples.Sensitivities import E2_fluxes3 as limits
from NuRadioMC.utilities import cross_sections

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from radiotools import plthelpers as php
from scipy import interpolate
from scipy import integrate
from scipy import interpolate as intp
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import coordinates as coord
from astropy.time import Time
from astropy.time import TimeDelta
import rnog_analysis_efficiency as rnogana
import pickle
import h5py
import glob
import yaml
import copy
import json
import os

n_deep = 12*12
n_shallow = 13*13 + n_deep

# BAC commented this out in Feb 2022 (since it's seemingly not used in the script)
# gen2_eng = []
# gen2_aeff = []
# f = open("gen2_optical_aeff.csv")
# for line in f:
#     if("#" in line):
#         continue
#     line_ = line.split(", ")
#     gen2_eng += [float(line_[0])]
#     gen2_aeff += [float(line_[1])]

# gen2_eng = np.power(10.0, np.array(gen2_eng).astype("float"))
# gen2_aeff = np.array(gen2_aeff).astype("float")
# gen2_veff = gen2_aeff * units.km * units.km * cross_sections.get_interaction_length(gen2_eng * units.eV, cross_section_type = 'ctw')
# gen2_veff_sr = gen2_veff * np.pi * 4.0

f = open("tabulated_veff_aeff_review_hybrid.csv")
eng = []
deep_veff = []
shallow_veff = []
overlap_fraction = []
deep_only = []
shallow_only = []

for line in f:
    if("#" in line):
        continue
    line_ = line.split(", ")
    if("energy" in line):
        for i in range(len(line_)):
            print(i, line_[i])
        continue

    eng += [float(line_[0])]
    deep_veff += [float(line_[1])]
    shallow_veff += [float(line_[3])]
    overlap_fraction += [float(line_[5])]
    deep_only += [float(line_[6])]
    shallow_only += [float(line_[7])]

eng = np.array(eng)
deep_veff = np.array(deep_veff)
shallow_veff = np.array(shallow_veff)
overlap_fraction = np.array(overlap_fraction)
deep_only = np.array(deep_only)
shallow_only = np.array(shallow_only)

fig, ax = limits.get_E2_limit_figure(show_ice_cube_EHE_limit=True,
                                     show_ice_cube_HESE_data=False,
                                     show_ice_cube_HESE_fit=False,
                                     show_ice_cube_mu=True,
                                     show_anita_I_III_limit=False,
                                     show_auger_limit=False,
                                     diffuse=True, 
                                     show_grand_10k=False, 
                                     show_grand_200k=False, 
                                     show_Heinze=False, 
                                     show_TA=False,
                                     show_ara=True, 
                                     show_arianna=True, 
                                     show_IceCubeGen2=False, 
                                     shower_Auger=False)
ax.grid(alpha=0.5)

livetime = 9.0 * units.year
total_veff = ((n_deep * deep_veff) + (n_shallow * shallow_veff))*(1/(1+overlap_fraction))
flux = fluxes.get_limit_e2_flux(energy=np.power(10.0, eng) * units.eV,
                                veff_sr= total_veff * units.km * units.km * units.km  * units.sr, 
                                livetime=livetime)
#ax.plot(np.power(10.0, eng) * units.eV / limits.plotUnitsEnergy, flux / limits.plotUnitsFlux, 'C0--')

veff_pess_deep = deep_veff * rnogana.PA_efficiency_pess(eng + 0.5)
veff_pess_shallow = shallow_veff * 0.79

veff_opp_deep = deep_veff * rnogana.PA_efficiency_opp(eng + 0.5)
veff_opp_shallow = shallow_veff * 1.0

total_veff_pess = ((n_deep * veff_pess_deep) + (n_shallow * veff_pess_shallow))*(1/(1+overlap_fraction))
total_veff_opp = ((n_deep * veff_opp_deep) + (n_shallow * veff_opp_shallow))*(1/(1+overlap_fraction))
flux_pess = fluxes.get_limit_e2_flux(energy=np.power(10.0, eng) * units.eV,
                                     veff_sr= total_veff_pess * units.km * units.km * units.km  * units.sr, 
                                     livetime=livetime)
flux_opp = fluxes.get_limit_e2_flux(energy=np.power(10.0, eng) * units.eV,
                                     veff_sr= total_veff_opp * units.km * units.km * units.km  * units.sr, 
                                     livetime=livetime)
SES = False
if(SES):
    ax.fill_between(np.power(10.0, eng) * units.eV / limits.plotUnitsEnergy, 
                    flux_pess / limits.plotUnitsFlux / 2.44, 
                    flux_opp / limits.plotUnitsFlux / 2.44, 
                    alpha=0.5,
                    color='red')

    legendfontsize = 11
    ax.annotate('IceCube Gen2 Radio SES', 
                xy=(1e19 * units.eV, 1.8e-10), xycoords='data',
                horizontalalignment='center', color='red', rotation=40, fontsize=legendfontsize)
    ax.set_ylim(1e-10, 1e-7)
    fig.savefig(f"DanSmith_projected_limit_SES.png", dpi=300)
else:
    ax.fill_between(np.power(10.0, eng) * units.eV / limits.plotUnitsEnergy, 
                    flux_pess / limits.plotUnitsFlux / 2.44, 
                    flux_opp / limits.plotUnitsFlux / 2.44, 
                    alpha=0.5,
                    color='red')

    legendfontsize = 11
    ax.annotate('IceCube Gen2 Radio 90% CL', 
                xy=(1e19 * units.eV, 1.5e-10), xycoords='data',
                horizontalalignment='center', color='red', rotation=40, fontsize=legendfontsize)

    ax.set_ylim(1e-10, 1e-7)
    fig.savefig(f"DanSmith_projected_limit_CL.png", dpi=300)

##################################################

fig, ax = limits.get_E2_limit_figure(show_ice_cube_EHE_limit=True,
                                     show_ice_cube_HESE_data=False,
                                     show_ice_cube_HESE_fit=False,
                                     show_ice_cube_mu=True,
                                     show_anita_I_III_limit=False,
                                     show_auger_limit=False,
                                     diffuse=True, 
                                     show_grand_10k=False, 
                                     show_grand_200k=False, 
                                     show_Heinze=False, 
                                     show_TA=False,
                                     show_ara=True, 
                                     show_arianna=True, 
                                     show_IceCubeGen2=False, 
                                     shower_Auger=False)
ax.grid(alpha=0.5)

livetime = 9.0 * units.year
total_veff = ((n_deep * deep_veff) + (n_shallow * shallow_veff))*(1/(1+overlap_fraction))
total_flux = fluxes.get_limit_e2_flux(energy=np.power(10.0, eng) * units.eV,
                                      veff_sr= total_veff * units.km * units.km * units.km  * units.sr, 
                                      livetime=livetime)

deep_veff = n_deep * deep_veff
deep_flux = fluxes.get_limit_e2_flux(energy=np.power(10.0, eng) * units.eV,
                                      veff_sr= deep_veff * units.km * units.km * units.km  * units.sr, 
                                      livetime=livetime)
shallow_veff = n_shallow * shallow_veff
shallow_flux = fluxes.get_limit_e2_flux(energy=np.power(10.0, eng) * units.eV,
                                        veff_sr= shallow_veff * units.km * units.km * units.km  * units.sr, 
                                        livetime=livetime)

SES = True
if(SES):
    ax.plot(np.power(10.0, eng) * units.eV / limits.plotUnitsEnergy, 
            total_flux / limits.plotUnitsFlux / 2.44,
            marker="D",
            alpha=1.0,
            color='red',
            label="Total Array")

    ax.plot(np.power(10.0, eng) * units.eV / limits.plotUnitsEnergy, 
            deep_flux / limits.plotUnitsFlux/ 2.44, 
            marker="o",
            alpha=1.0,
            color='C0',
            label="Deep Component, 144 Stations")

    ax.plot(np.power(10.0, eng) * units.eV / limits.plotUnitsEnergy, 
            shallow_flux / limits.plotUnitsFlux/ 2.44, 
            marker="D",
            alpha=1.0,
            color='C1',
            label="Shallow Component, 313 Stations")

    ax.legend(loc="upper right")

    legendfontsize = 11
    '''
    ax.annotate('IceCube Gen2 Radio SES', 
                xy=(1e19 * units.eV, (1.8*2.0)*1e-10), xycoords='data',
                horizontalalignment='center', color='red', rotation=40, fontsize=legendfontsize)
    '''
    ax.annotate('IceCube Gen2 Radio SES', 
                xy=(2.7e16 * units.eV, (1.5)*1e-10), xycoords='data',
                horizontalalignment='center', color='red', rotation=-65, fontsize=legendfontsize)
    ax.set_ylim(1e-10, 1e-7)
    fig.savefig(f"DanSmith_projected_limit_SES_trigger_level.png", dpi=300)
else:
    ax.plot(np.power(10.0, eng) * units.eV / limits.plotUnitsEnergy, 
            total_flux / limits.plotUnitsFlux,
            marker="D",
            alpha=1.0,
            color='red',
            label="Total Array")

    ax.plot(np.power(10.0, eng) * units.eV / limits.plotUnitsEnergy, 
            deep_flux / limits.plotUnitsFlux,
            marker="o",
            alpha=1.0,
            color='C0',
            label="Deep Component, 144 Stations")

    ax.plot(np.power(10.0, eng) * units.eV / limits.plotUnitsEnergy, 
            shallow_flux / limits.plotUnitsFlux,
            marker="D",
            alpha=1.0,
            color='C1',
            label="Shallow Component, 313 Stations")

    ax.legend(loc="lower left")
    #legendfontsize = 11
    ax.annotate('90% CL', 
                xy=(1e19 * units.eV, 2.8e-10), xycoords='data',
                horizontalalignment='center', color='red', rotation=40, fontsize=legendfontsize)

    ax.set_ylim(1e-10, 1e-7)
    fig.savefig(f"DanSmith_projected_limit_CL_trigger_level.png", dpi=300)
