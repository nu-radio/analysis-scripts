# combine the standard detectors
CHANNEL_FILE="single_station/Gen2_trigger_channels.json"
for STATION_FILE in "Gen2_baseline_array.json" "Gen2_hex_shallow_array.json" "Gen2_hex_shallowheavy_array.json" "Gen2_hex_hybrid_only_array.json" "Gen2_hex_shallow_1.2km_array.json"
do
    python3 scripts/combine_channels_stations.py --output="trigger_$STATION_FILE" --stations_from="array_layouts/$STATION_FILE" --channels_from="$CHANNEL_FILE"
done

# combine the HPol trigger arrays
CHANNEL_FILE="single_station/Gen2_4channel_HPol_VPol_trigger_channels.json"
for STATION_FILE in "Gen2_baseline_array.json" "Gen2_hex_shallow_array.json" "Gen2_hex_shallowheavy_array.json" "Gen2_hex_hybrid_only_array.json" "Gen2_hex_shallow_1.2km_array.json"
do
    python3 scripts/combine_channels_stations.py --output="hpol_trigger_$STATION_FILE" --stations_from="array_layouts/$STATION_FILE" --channels_from="$CHANNEL_FILE"
done
