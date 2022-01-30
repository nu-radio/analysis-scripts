# Gen2 Radio Simulations
repository for gen2 radio sim production stuff

## Overview

The strategy for the simulation requires a few stages:

- step 0: create the configuration files for how to distribute the vertices, zeniths, azimuths, energies, etc of incident neutrinos
- step 1: generate a list of energy depositions in the ice based on the neutrino primaries from step 0 (from neutrino primary interactions, but also secondary interactions like pair production, etc. in the case of electron and muon neutrinos)
- step 2: propagate radio emission from the energy depositions in step 1 to a simulated station, run triggering, etc.
- step 3: combine all of the small "output" files into a single main file

The simulation inputs, Step 0 and 1, are defined by input configuration names such as `secondaries_1700km2`, which suggest that principle fiducial volume is distributed over a 1700 km^2 area. For the tau and muon simulation, the volume is automatically enlarged by NuRadioMC.

The simulation outputs, Step 2 and 3, are defined by a combination of:
1. a geometry specification (how antennas are distributed in the ice, etc.), defined in a `.json` file
2. a configuration specification (what Askaryan and ice model to use, etc)  defined in a `.yaml` file
3. a simulation configuration (what triggers to run, etc.) defined in a `.py` file.
The outputs need to be organized according to these three specification files.

Because the inputs (step 0 and 1) are totally separable from the outputs (step 2 and 3), they should be stored in a separate directory. The structure used is below:

```
simulation_input
     secondaries_1700km2
          step0
               e
               mu
               tau
          step1
               e
               mu
               tau
simulation_output
     secondaries_1700km2
          geometry_specification
               configuration_specification
                    simulation_specification
                         e
                         mu
                         tau
```

We are also running the simulation in specific zenith bands, and so for each `e`, `mu`, and `tau` folder, there is a zenith band folder, e.g.: `e_18.50eV_-0.2_0.1` for electron neutrinos, of primary energy 18.50eV, and zeniths drawn from the range -0.2 to -0.1.

Simulations are generated in "parts"--that is, for a given flavor, energy bin, and zenith bin, there might be 10 "parts" which must be run in parallel, and then merged in step 3.

This directory for production scripts only includes the scripts for step 0, which are "cluster agnostic." Step 1-3 will depend on whether the user is using an htcondor cluster, or slurm cluster, etc. are are not included in this repository.

## Simulation Tools and Steps

We are using [NuRadioMC tag v2.1.3](https://github.com/nu-radio/NuRadioMC/releases/tag/v2.1.3) and tag [NuRadioMC tag v2.1.5](https://github.com/nu-radio/NuRadioMC/releases/tag/v2.1.5)
The >=X.X.3 hotfix version is important, as during TDR simulation prep, a bug was found in the units used in the proposal shower energy cut.
(See [this](https://github.com/nu-radio/NuRadioMC/pull/363) pull request.)

More explanation on why multiple versions were used is explained below.

There is a small "helper" class provided in this repository. This standardizes things like the zenith binning, energy binning, etc. To use it, you will need to do:

`export PYTHONPATH=/path/to/tools/directory/:$PYTHONPATH`

### NuRadioMC Versions
We realized mid-production that there was a bug in versions less than v2.1.5. 
The bug is that the trigger time was not set correctly in the phased array
module, and so the time traces were cut incorreclty in the .nur files.
(Meaning the waveforms are trimmed incorrectly and 
the pulse is often not present in the simulated waveform.)
The hdf5 files, and the veffs, coincidences, etc, are all correct.
This bug is discussed more [here](https://github.com/nu-radio/NuRadioMC/pull/372).

Because we distributed production, and did the arrays out of order, different
versions were used on differente arrays and energy ranges.
Again, this should only affect .nur files, not the .hdf5 files 
and the general high-level physics conclusions.

Information on which versions were used where is below. 

**1E19 and Above (computed by Brian on the Grid)**
- Baseline, Hybrid-Only:  v2.1.4
- Shallow, Shallow-Heavy: v2.1.5

**1E18.5 and Below (computed by Steffen on the Desy cluster)**
- Baseline, Hybrid-Only, Shallow, Shallow-Heavy:  v2.1.4


This information is also contained in the hdf5 files themselves,
and so can be recovered after the fact.  It can be extracted by doing:
```python
f = h5py.open('filename.hdf5')
f.attrs['NuRadioMC_version']
```


## step 0 scripts
Step 0 is where we output the "control" `.py` files for NuRadioMC that tell NuRadioMC how to distribute the neutrinos that will become the step 1 files. This includes the neutrino energy, zeniths, vertex volume, etc. Run this by doing:

`python make_step0_files.py`

The main parameter, which sits near the top of the file, is the `base_dir`, which is where you ant the step 0 and step 1 files to be placed.
