# Readme

## Radio Files

The radio files contained here are from the Feb 2021 review array.

`review_array_dict_e.pkl` contains the radio effective volumes for the array,
in a style expected by by NuRadioMC.
Same with `mu` and `tau`.

The `tabulated_veff_aeff_review_hybrid.csv` are the array volumes
broken out by surface and deep. Using this files requires applying scaling factors.

## Optical Files

The optical files are provided by Max Meier (Chiba U).
Max devied a pseudo-event selection based on Gen2-Optical simulations.
So these are **analysis level** numbers.

From Max: 

The files Gen2_EHE_effective_area_{e|mu|tau}.npz hold the effective area
for a simple EHE analysis port to IceCube Gen2 using simulations with the
Sunflower 240m detector configuration and pDOMs as detector modules (essentially
IceCube DOMs with HQE PMTs).

The nu/nubar averaged effective area, given in m2, is binned in three dimensions:
    - True cos(zenith)
    - Energy at the detector volume
    - Neutrino energy at the surface

The binning scheme used is given in the npz file.

In short, to the effective area vs energy and zenith, a user should sum over
the second dimension. To the get the sky-averaged effective area,
the user should average over the first dimension, and multiply by pi.
This would produce aeff*str.