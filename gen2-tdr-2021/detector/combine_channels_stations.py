import numpy as np
import json
import argparse
"""
a script to combine a json dict containing only the stations with the json file where channels are defined
"""
parser = argparse.ArgumentParser()
parser.add_argument('--stations_from', type=str, required=True, help="json file to pick stations from")
parser.add_argument('--channels_from', type=str, required=True, help="json file to pick channels from")
parser.add_argument('--output', type=str, default="merged_channels_stations.json")

args = parser.parse_args()
print(args)
with open(args.channels_from) as f:
    channel_data = json.load(f)
with open(args.stations_from) as f:
    station_data = json.load(f)

channel_data["stations"] = station_data

with open(args.output, 'w') as json_file:
  json.dump(channel_data, json_file, indent=4)

