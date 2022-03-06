
# you need the plot_helper (from one directory up) in your python path
# you also need to set the NOTEBOOKDIR (for one directory up)
# this needs to be an absolute path, sorry...
export PYTHONPATH=$PWD/..:$PYTHONPATH
export NOTEBOOKDIR="/home/brianclark/Gen2/radio/analysis-scripts/gen2-tdr-2021/analysis_notebooks"

declare -a dets=("0" "1" "2" "3")

#############################
# first, the deep/shallow and hybrid/shallow mode calculation
#############################

declare -a modes=("deepshallow" "hybridshallow")
for det in "${dets[@]}"
do
    for mode in "${modes[@]}"
    do
        python calc_fractions.py --det $det --mode $mode # first, calculate the information
        python get_veffs.py --det $det --mode $mode # now, make plots and write to csv
    done
done


#############################
# now we need to get the review array
#############################

python plot_veffs.py --det 4 --mode "review"


#############################
# and now to make plots
#############################

python plot_stats.py
python plot_limits.py


#############################
# now, the advanced mode statistics calculations
#############################

declare -a modes=("stations_and_triggers", "triggers_coincidences")
for det in "${dets[@]}"
do
    for mode in "${modes[@]}"
    do
        python calc_fractions_advanced.py --det $det --mode $mode # first, calculate the information
        python get_veffs_advanced.py --det $det --mode $mode # now, make plots and write to csv
    done
done



