import numpy as np
from NuRadioReco.utilities import units
import helper as hp

flavors = [
    "e", 
    "mu", 
    "tau"
    ]

coszenbins = hp.get_coszenbins()
logEs = hp.get_logEs()
logEs = np.arange(16, 20.1, 0.5)
energies = 10 ** logEs * units.eV

step3dir = "/data/sim/Gen2/radio/2020/gen2-tdr-2021/simulation_output/secondaries_1700km2/step3/"
step4dir = "/data/sim/Gen2/radio/2020/gen2-tdr-2021/simulation_output/secondaries_1700km2/step4/unmerged/"

det_files_dict = {
    "gen2r_baseline" : "baseline_array",
    "gen2r_hybrid" : "hex_hybrid_only_array",
    "gen2r_shallow" : "hex_shallow_array",
    "gen2r_shallowheavy" : "hex_shallowheavy_array",
}
det_files_labels = [
    "gen2r_baseline", 
    "gen2r_hybrid" , 
    "gen2r_shallow", 
    "gen2r_shallowheavy"
    ]
config_file = "config_ARZ2020_noise"
sim_file = "D01detector_sim"

dag_file_name='dagman_step4.dag'
instructions = ""
instructions += 'CONFIG config.dagman\n\n'
with open(dag_file_name, 'w') as f:
    f.write(instructions)

master_index=0

for det_file_label in det_files_labels:

    det_file = det_files_dict[det_file_label]

    for flavor in flavors:

        for iE in range(len(logEs)):

            for iC in range(len(coszenbins)-1):
                
                czen1 = coszenbins[iC]
                czen2 = coszenbins[iC + 1]
                E = energies[iE]

                instructions = ""
                instructions += f'JOB job_{master_index} step4_job_standardize.sub \n'
                instructions += f'VARS job_{master_index} step3dir="{step3dir}" step4dir="{step4dir}" detfile="{det_file}" configfile="{config_file}" simfile="{sim_file}" flavor="{flavor}" energy="{logEs[iE]:.2f}" czmin="{czen1:.1f}" czmax="{czen2:.1f}" \n\n'

                with open(dag_file_name, 'a') as f:
                    f.write(instructions)

                master_index+=1


