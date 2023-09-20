## Step 4

### Introduction

These are the scripts for making the step 4 files.
This step is meant to be run on the UW-Madison cluster
on `submit` (not `sub-1`, the grid machine).
So you need access to the local filesystem for this to work.

What happens at step 4 is all the merged hdf5 files
are converted into json files, which are easier to read.
These create one json file for each energy and zenith bin.
These json files are then merged.

Note: the code "cleverly" creates json files full of zeros 
(with the `dummy.json` file as a template), so that the result of the merging
always has the right number of energy and zenith bins, regardless of
if a file was missing or empty.


### How to Use

Start by making the dagman file:

`python make_step4_dag.py`

Then, you can just submit them:

`condor_submit_dag dagman_step4.dag`

Be careful to change any of the paths (that are probably relative
to Brian at the moment).

If you want to change if things are getting over-sampled or not,
go into `step4.py` and change things like:
`make_dummy` to `make_dummy_oversampling`.
Also, don't forget to change the `oversampling_theta` argument in 
`get_Veff_Aeff`.