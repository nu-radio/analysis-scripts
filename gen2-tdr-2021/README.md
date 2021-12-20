# config
Contains a config file for the simulation. Differences to February review: Uses ARZ2020 for raytracing and PREM for weight calculation.

# detector
Contains sub-directories
* `antenna_arrangement`: holding antenna layounts for only one `shallow-only` and `hybrid` stations. I.e. positions of all (triggering) antennas
* `station_positions`: the array layouts to be simulated. Only the stations (positions, reference stations) are defined there). Reference stations are `1001` (hybrid) and `2001` (shallow-only)
The script `combine_channels_stations` can be used to combine antenna layouts with station positions. cf. the `.sh` file for details.

# detsim
Contains the python script for trigger simulation.
