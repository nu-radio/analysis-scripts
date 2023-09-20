import os
import h5py
import numpy as np
from glob import glob
from NuRadioMC.utilities import merge_hdf5
from NuRadioMC.utilities import Veff
import sys
from helper import get_nsim_ntrig
from helper import trigger_combinations

""" script that calculates Veff json files while merging on the fly, adds a set of combined (OR) triggers for deep/shallow  """

import argparse

parser = argparse.ArgumentParser(description='Generate effective volumes **per non-merged directory** with additional stats on simulated/triggered events (doing merging on the fly)')

parser.add_argument('temp_dir', type=str, 
                    help='temporary (step3) directory to store merged files and delete them afterwards')
parser.add_argument('--pattern', type=str, default="/",
                    help='pattern in the directory path to match, necessary if running multiple simultaneous jobs on the step2 dir.') #default should always match

parser.add_argument('--step2dir', type=str, default="/lustre/fs22/group/radio/shallman/gen2-tdr-2021/simulation_output/secondaries_1700km2/step2",
                    help='step2 directory (input directory) of simulations to walk through')
parser.add_argument('--step4dir', type=str, default="/lustre/fs22/group/radio/shallman/gen2-tdr-2021/simulation_output/secondaries_1700km2/step4",
                    help='step4 directory (output directory) to store per-directory Veff files')

args = parser.parse_args()


step2dir = args.step2dir
step3dir = args.temp_dir
step4dir = args.step4dir
patternmatch = args.pattern

print(f"adding combined triggers to veffs: {trigger_combinations}")

step2_dirs = glob(f"{step2dir}/**/", recursive = True)
print(f"found {len(step2_dirs)} step2 dirs in {step2dir}")


for directory in step2_dirs:
    if patternmatch is not None:
        if not patternmatch in directory:
            continue

    h5files = glob(f'{directory}*.hdf5')
    if not h5files:
        print("no hdf5 files in path, continuing")
        continue
    else:
        # check if json already exists
        identifier = directory.replace(step2dir, "")
        identifier = "__".join([x for x in identifier.split("/") if x != ""])
        output_jsonfile = f"Veff__{identifier}.json"
        if os.path.exists(os.path.join(step4dir, output_jsonfile)):
            print(f"Veff json file already exists")
            continue



        print(directory)
        outdir = directory.replace(step2dir, step3dir)
        while outdir.endswith("/"):
            outdir = outdir[:-1]
        outfile = outdir.split("/")[-1] + ".hdf5"
        outdir = os.path.abspath(os.path.join(outdir, os.pardir))
        if not os.path.exists(step3dir):
            os.makedirs(step3dir)
        
        h5files.sort()
        print(f"- OUTDIR: {outdir}")
        print(f"- OUTFILE: {outfile}")

        if os.path.exists(os.path.join(step3dir, outfile)):
            print(f"merged file already exists")
            continue
        temp_out = os.path.join(step3dir, outfile)
        merge_hdf5.merge2(h5files, temp_out)
        # now get the stats

    if not os.path.isfile(temp_out):
        print("no hdf5 files in path, continuing")
        continue
    else:

        #identifier = directory.replace(step2dir, "")
        #identifier = "__".join([x for x in identifier.split("/") if x != ""])
        #output_jsonfile = f"Veff__{identifier}.json"

        if not os.path.exists(step4dir):
            os.makedirs(step4dir)

        print(f"- OUTDIR: {step4dir}")
        print(f"- OUTFILE: {output_jsonfile}")

        nsim, ntrig = get_nsim_ntrig(temp_out)
        print("SIM/TRIG", nsim, ntrig)
        if ntrig == 0:
            print(f"no events triggered")
            os.remove(temp_out) 
            continue
        v = Veff.get_Veff_Aeff(step3dir, trigger_combinations=trigger_combinations)
        if not len(v)==1:
            print("something seems to have gone wrong. length not one")
    
        v[0]["stats_ntrig"] = ntrig
        v[0]["stats_nsim"] = nsim
        v[0]["nfiles"] = len(h5files)
        Veff.export(os.path.join(step4dir, output_jsonfile), v, export_format="json")
    os.remove(temp_out)    
    #break
