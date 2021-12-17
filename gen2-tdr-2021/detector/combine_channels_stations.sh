CHANNEL_FILE="antenna_arrangement/Gen2_shallow_and_PA_only.json"
for STATION_FILE in "Gen2_baseline_array.json" "Gen2_hex_shallow_array.json" "Gen2_hex_shallowheavy_array.json" "Gen2_hybrid_only_array.json"
do
    python combine_channels_stations.py --output="trigger_$STATION_FILE" --stations_from="station_positions/$STATION_FILE" --channels_from="$CHANNEL_FILE"
done
