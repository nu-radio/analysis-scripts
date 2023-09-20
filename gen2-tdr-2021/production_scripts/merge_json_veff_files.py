from glob import glob
import pandas as pd
import json
import sys
import os

import argparse

parser = argparse.ArgumentParser(description='Join together individual Veff files per E and cosz to one file.')
parser.add_argument('--input_dir', type=str, default = ".",
                    help='input directory')
parser.add_argument('--output_dir', type=str, default = ".",
                    help='output directory')
parser.add_argument('--separator', type=str, default = "__",
                      help='output directory')


args = parser.parse_args()

indir = args.input_dir
outdir = args.output_dir

#if len(sys.argv)>1:
#    indir = sys.argv[1]
#if len(sys.argv)>2:
#    outdir = sys.argv[2]    

files = glob("/".join([indir, "*.json"]))

def get_identifiers(files, separator):
    identifiers = []
    for f in files:
        base = os.path.basename(f)
        if not "eV_" in base:
            # already a merged file
            continue
        identifiers.append(base[:base.rfind(separator)])
    identifiers = list(set(identifiers))
    return identifiers

def matching_files(files, identifier):
    matches = []
    for f in files:
        if identifier in f:
            matches.append(f)
    return matches 

identifiers = get_identifiers(files, args.separator)
print("identifiers are {}".format(identifiers))

for identifier in identifiers:
    mfiles = matching_files(files, identifier)
    if len(mfiles) == 0:
        continue
    print("concatenating", identifier, len(mfiles))
    #for m in mfiles:
    #    print(f". {m}")
    data = [pd.read_json(f) for f in mfiles]
    cdata = pd.concat(data, ignore_index=True)
    cdata.sort_values(["energy", "thetamax"], ascending=[True, False], ignore_index=True)
    out = "/".join([outdir, identifier+".json"])
    print(f"-> OUTPUT: {out}")
    cdata.to_json(out, orient="records", indent=4)

