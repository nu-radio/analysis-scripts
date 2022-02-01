# Readme

The code in this directory make plots of plots.

It's main goal is to calculate the veff and aeff,
but broken out in a more sophistocated fraction.
Meaning deep vs shallow trigger, hybrid vs surface stations, 
and all the different trigger components (hybrid station, deep trigger vs 
shallow station surface trigger, etc).
See calc_fractions.py and calc_fractions_advanced.py 
for more discussion.

It can also make plots of the triggering statistics,
plots of the veff/aeff, and projected limits.

To get all the plots at once, run `bash do_calculations.sh`.

A fair warning, it takes about 20-30 minutes, as it has to 
go over all the hdf5 files for all arrays (4), all energies (9),
and all zenith bins (20). It is multithreaded, so reasonably fast.

The resuls of the triggering calculation are stored in the `data` folder.
The plots are dumped in `plots`, and the tabulated veff/aeff
are put into `results`.