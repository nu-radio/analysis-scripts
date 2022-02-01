#!/bin/bash

# load variables
step3dir=$1
step4dir=$2
detfile=$3
configfile=$4
simfile=$5
flavor=$6
energy=$7
czmin=$8
czmax=$9

export starting=$PWD

merged_file=${flavor}_${energy}eV_${czmin}_${czmax}.hdf5
out_file='Veff__'${detfile}'__config_'${configfile}'__'${simfile}'__'${flavor}'__'${energy}'eV_'${czmin}'_'${czmax}'.json'

# in production
working_dir=$TMPDIR

# clear out the python path and source the setup file for gen2 radio simulations
unset PYTHONPATH
ls /cvmfs/icecube.opensciencegrid.org/users/brianclark/gen2radiosim/setup.sh
source /cvmfs/icecube.opensciencegrid.org/users/brianclark/gen2radiosim/setup.sh

# need the helper
export PYTHONPATH=/home/brianclark/Gen2/radio/analysis-scripts/gen2-tdr-2021/production_scripts:$PYTHONPATH

export HDF5_USE_FILE_LOCKING='FALSE'
cp /home/brianclark/Gen2/radio/gen2_radio_sims/scripts/gen2-tdr-2021/step4/step4_standardize.py .
cp /home/brianclark/Gen2/radio/gen2_radio_sims/scripts/gen2-tdr-2021/step4/dummy.json .
python step4_standardize.py ${step3dir}/${detfile}/${configfile}/${simfile}/${flavor}/${merged_file} ${out_file} ${energy} ${czmin} ${czmax}
rm dummy.json

if test -f "$out_file"; then
    echo "The output json file -- $out_file -- was made correctly. Continue with transferring"

    mv *.json ${step4dir}/.
else
    echo "The output json file -- $out_file -- was NOT made correctly. Report failure"
    exit 1
fi

cd $starting
