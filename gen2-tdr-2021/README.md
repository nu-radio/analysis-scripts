# config
Contains a config file for the simulation. Differences to February review:
* ARZ2020 for raytracing
* PREM for weight calculation (three layer model has rock up to surface, which might change weights a bit around the horizon)
* events are split up if time differences are larger than 1000ns (instead of 100ns). This is needed for the simulation of a hybrid station with time delays between the shallow and deep components. 

# detector
Contains the 4 json detector description files for the 4 arrays under study. 

There are two bash scripts, which can be used to (re)build the `array_layouts` for the arrays to be simulated and to combine the trigger channels with the array layouts to have `.json` files for running the trigger.

Contains sub-directories
* `scripts`: directory containing python scripts to generate `array_layouts` and to combine the channels from `single_station` with `array_layouts`
* `single_station`: holding antenna layouts for only one `shallow-only` and `hybrid` station each. I.e. positions of all (triggering) antennas
* `array_layouts`: the array layouts to be simulated. Only the stations (positions, reference stations) are defined there). Reference stations are `1001` (hybrid) and `2001` (shallow-only). 

# detsim
Contains the python script for trigger simulation.
There is only one script that simulated both the LPDA surface trigger as well as the 4 and 8 channel deep trigger. Based on the number of stations it is automatically determined if it is a hybrid station or shallow only station. 
The trigger channels are low pass filtered on-the-fly to the lower trigger bandwidth. The resulting noise RMS is recalculated so that the trigger thresholds are correctly set. 
