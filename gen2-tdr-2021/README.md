# config
Contains a config file for the simulation. Differences to February review:
* ARZ2020 for raytracing
* PREM for weight calculation (three layer model has rock up to surface, which might change weights a bit around the horizon)

# detector
There are two bash scripts, which can be used to build the `array_layouts` for the arrays to be simulated and to combine the trigger channels with the array layouts to have `.json` files for running the trigger.


Contains sub-directories
* `scripts`: directory containing python scripts to generate `array_layouts` and to combine the channels from `single_station` with `array_layouts`
* `single_station`: holding antenna layouts for only one `shallow-only` and `hybrid` station each. I.e. positions of all (triggering) antennas
* `array_layouts`: the array layouts to be simulated. Only the stations (positions, reference stations) are defined there). Reference stations are `1001` (hybrid) and `2001` (shallow-only). 

# detsim
Contains the python script for trigger simulation.
