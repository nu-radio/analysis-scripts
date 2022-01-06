import numpy as np
import os

from NuRadioReco.utilities import units
import helper as hp

# this is the "base" directory for TDR simulations
base_dir = "/data/sim/Gen2/radio/2020/gen2-tdr-2021/simulation_input/"

# make necessary step 0 and step 1 subdirectories (if they don't already exist)
step0dir = os.path.join(base_dir, f"secondaries_1700km2", "step0")
if(not os.path.exists(step0dir)):
    os.makedirs(step0dir)

step1dir = os.path.join(base_dir, f"secondaries_1700km2", "step1")
if(not os.path.exists(step1dir)):
    os.makedirs(step1dir)

coszenbins = hp.get_coszenbins()
logEs = hp.get_logEs()
energies = 10 ** logEs * units.eV

phimin = 0.*units.deg
phimax = 360.*units.deg

# these are the fiducial volume settings that give ~5km clearange
# on from the edge of the largest array under consideration 
# (which is the the shallow heavy design)
volume = {'fiducial_xmin': -16. * units.km,
        'fiducial_xmax': 18. * units.km,
        'fiducial_ymin': -25. * units.km,
        'fiducial_ymax': 25. * units.km,
        'fiducial_zmin': -2.7 * units.km,
        'fiducial_zmax': 0,
        'full_zmin': -3.3 * units.km,
        'full_zmax': 0
        }

distance_cut_polynomial = np.polynomial.polynomial.Polynomial([-1.56610502e+02, 2.54131322e+01, -1.34932379e+00, 2.39984185e-02])

def get_distance_cut(shower_energy):
    return max(100 * units.m, 10 ** distance_cut_polynomial(np.log10(shower_energy)))

flavors = [
            "e", 
            "mu", 
            "tau"
            ]

flavor_ids = {
            'e': [12,-12],
            'mu': [14,-14],
            'tau': [16,-16]
            }

for flavor in flavors:

    run_proposal = 'True'
    if flavor == "e":
        # if it's nu-e events, we can turn proposal off
        run_proposal = 'False'

    for iE in range(len(logEs)):


        max_dist = get_distance_cut(10 ** logEs[iE])
        print(f"maximum radius for E = {10**logEs[iE]:.2g} is {max_dist/units.m:.0f}m")

        # the ice sheet at SP is 2.7km deep, we add 200meters to be save for 200m deep dipole simulations
        volume['fiducial_zmin'] = -min(max_dist + 200*units.m, 2.7 * units.km)  
        
        for iC in range(len(coszenbins) - 1):
            czen1 = coszenbins[iC]
            czen2 = coszenbins[iC + 1]
            E = energies[iE]
            minEproposal = 1E-3 * E
            thetamax = np.arccos(czen1)
            thetamin = np.arccos(czen2)
            pattern = f"{flavor}_{logEs[iE]:.2f}eV_{czen1:.1f}_{czen2:.1f}"
            num_parts, num_events = hp.get_number_of_parts_and_events(flavor, logEs[iE], czen1)
            print(pattern)
            
            folder = os.path.join(step0dir, flavor, f"{pattern}")
            if(not os.path.exists(folder)):
                os.makedirs(folder)

            folder1 = os.path.join(step1dir, flavor, f"{pattern}")
            if(not os.path.exists(folder1)):
                os.makedirs(folder1)

            for ijob in range(num_parts):
                
                instruction = ""
                instruction += 'from NuRadioMC.EvtGen.generator import generate_eventlist_cylinder\n'
                instruction += 'from NuRadioReco.utilities import units\n\n'
                out_filename = "in_" + f"{pattern}" + f".part{ijob:06}" + ".hdf5"
                start_event_id = int(ijob * (num_events*5)) + 1
                instruction += f"generate_eventlist_cylinder('{out_filename}', {num_events}, {E}, {E}, {volume},\n"
                instruction += f"thetamin={thetamin}, thetamax={thetamax}, phimin={phimin},\n"
                instruction += f"phimax={phimax}, \n"
                instruction += f"start_event_id={start_event_id},\n"
                instruction += f"proposal={run_proposal}, proposal_config='config_PROPOSAL.json', n_events_per_file=None,\n"
                instruction += f"flavor={flavor_ids[flavor]},\n"
                instruction += f"proposal_kwargs={{'low_nu': {minEproposal}*units.eV, 'min_energy_loss_nu': {minEproposal}*units.eV}})\n"
                instruction += "\n"
                python_filename = f'{pattern}_{ijob:06d}.py'
                
                if(not os.path.exists(os.path.join(step0dir, flavor, pattern))):
                    os.makedirs(os.path.join(step0dir, flavor, pattern))

                if(not os.path.exists(os.path.join(step1dir, flavor, pattern))):
                    os.makedirs(os.path.join(step1dir, flavor, pattern))
                
                with open(os.path.join(step0dir, flavor, pattern, python_filename), 'w') as f:
                    f.write(instruction)
