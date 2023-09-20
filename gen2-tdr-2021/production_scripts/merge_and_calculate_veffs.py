import os
import numpy as np
from glob import glob
from NuRadioMC.utilities import merge_hdf5
from NuRadioMC.utilities import Veff
#from helper import get_nsim_ntrig
from helper import trigger_combinations

step2dir="/lustre/fs22/group/radio/shallman/gen2-tdr-2021/simulation_output/secondaries_1700km2/step2"
step3dir="/lustre/fs22/group/radio/shallman/gen2-tdr-2021/simulation_output/secondaries_1700km2/step3"
step4dir="/lustre/fs22/group/radio/shallman/gen2-tdr-2021/simulation_output/secondaries_1700km2/step4"

do_merging = False
do_veffs = True

if do_merging:
    step2_dirs = glob(f"{step2dir}/**/", recursive = True)
    print(f"found {len(step2_dirs)} step2 dirs in {step2dir}")

    for directory in step2_dirs:
        h5files = glob(f'{directory}*.hdf5')
        if not h5files:
            print("no hdf5 files in path, continuing")
            continue
        else:
            print(directory)
            outdir = directory.replace(step2dir, step3dir)
            while outdir.endswith("/"):
                outdir = outdir[:-1]
            outfile = outdir.split("/")[-1] + ".hdf5"
            outdir = os.path.abspath(os.path.join(outdir, os.pardir))
            if not os.path.exists(outdir):
                os.makedirs(outdir)
        
            h5files.sort()
            print(f"- OUTDIR: {outdir}")
            print(f"- OUTFILE: {outfile}")

            if os.path.exists(os.path.join(outdir, outfile)):
                print(f"merged file already exists")
                continue
            merge_hdf5.merge2(h5files, os.path.join(outdir, outfile))
            #break


if do_veffs:
    step3_dirs = glob(f"{step3dir}/**/", recursive = True)
    print(f"found {len(step3_dirs)} step3 dirs in {step3dir}")

    for directory in step3_dirs:
        h5files = glob(f'{directory}*.hdf5')
        if not h5files:
            print("no hdf5 files in path, continuing")
            continue
        else:
            print(directory)

            identifier = directory.replace(step3dir, "")
            identifier = "__".join([x for x in identifier.split("/") if x != ""])
            output_jsonfile = f"Veff__{identifier}.json"

            if not os.path.exists(step4dir):
                os.makedirs(step4dir)
    
            print(f"- OUTDIR: {step4dir}")
            print(f"- OUTFILE: {output_jsonfile}")

            if os.path.exists(os.path.join(step4dir, output_jsonfile)):
                print(f"Veff json file already exists")
                continue

            v = Veff.get_Veff_Aeff(step3dir, trigger_combinations=trigger_combinations)

            Veff.export(os.path.join(step4dir, output_jsonfile), v, export_format="json")
            #break


