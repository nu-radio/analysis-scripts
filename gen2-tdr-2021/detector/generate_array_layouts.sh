for layout in "baseline" "hex_hybrid_only" "hex_shallow" "hex_shallowheavy" "hex_shallow_1.2km"
do
    python3 "scripts/generate_${layout}_array_json.py" && mv "${layout}_array.json" "array_layouts/Gen2_${layout}_array.json"
done
