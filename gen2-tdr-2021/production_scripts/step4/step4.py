import numpy as np
import sys
import pickle as pickle
import os
import json
from helper import trigger_combinations
from NuRadioMC.utilities import Veff
import copy

def make_dummy(outfile, energy, czmin, czmax):
    thetamin = np.arccos(float(czmax)) # !! thetamin corresponds to czmax
    thetamax = np.arccos(float(czmin))
    energy = float(energy)

    # open the dummy file full of zeros
    f = open('dummy.json')
    data = json.load(f)
    the_data = data[0]

    # swap out the relevant keys
    the_data['thetamin'] = thetamin
    the_data['thetamax'] = thetamax
    the_data['energy'] = energy
    the_data['energy_max'] = energy
    the_data['energy_min'] = energy
        
    # write the result out
    with open(outfile, 'w') as outfile:
        json.dump([the_data], outfile, sort_keys=True, indent=4)

def make_dummy_oversampling(outfile, energy, czmin, czmax, oversampling):
    thetamin = np.arccos(float(czmax)) # !! thetamin corresponds to czmax
    thetamax = np.arccos(float(czmin))
    energy = float(energy)
    costhetamin = np.cos(thetamin)
    costhetamax = np.cos(thetamax)
    thetas = np.arccos(np.linspace(costhetamin, costhetamax, oversampling + 1))
    thetas_min = thetas[:-1]
    thetas_max = thetas[1:]

    # open the dummy file full of zeros
    f = open('dummy.json')
    data = json.load(f)
    the_data = data[0]

    to_write = []

    for i, t_min in enumerate(thetas_min):
        the_data_copy = copy.deepcopy(the_data)
        the_data_copy['thetamin'] = t_min
        the_data_copy['thetamax'] = thetas_max[i]
        the_data_copy['energy'] = energy
        the_data_copy['energy_max'] = energy
        the_data_copy['energy_min'] = energy
        the_data_copy['domega'] = the_data_copy['domega']/oversampling
        to_write.append(the_data_copy)
    
    # write the result out
    with open(outfile, 'w') as outfile:
        json.dump(to_write, outfile, sort_keys=True, indent=4)


if __name__ == "__main__":

    infile = sys.argv[1]
    outfile = sys.argv[2]
    logenergy = float(sys.argv[3])
    energy = 10 ** logenergy
    czmin = float(sys.argv[4])
    czmax = float(sys.argv[5])

    if not os.path.exists(infile):
        print("input file does not exist {}. Need to make a dummy file!".format(infile))
        # make_dummy(outfile, energy, czmin, czmax)
        make_dummy_oversampling(outfile, energy, czmin, czmax, 4)
        sys.exit()


    try:
        v = Veff.get_Veff_Aeff(infile, trigger_combinations=trigger_combinations,
            oversampling_theta=4
            )
        Veff.export(outfile, v, export_format="json")
    except:
        print("our attempt at actually making the file failed. Make a dummy.")
        # make_dummy(outfile, energy, czmin, czmax)
        make_dummy_oversampling(outfile, energy, czmin, czmax, 4)

