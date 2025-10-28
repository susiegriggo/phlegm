#!/usr/bin/env bash
set -euo pipefail

# Number of subunits to generate
max_subunits=10

# Input directory containing .a3m files
input_search="$1"

# Output directory
output_dir="$2"

# Create output directory
mkdir -p "$output_dir"

# Loop through each .a3m file in the input directory
for file in "$input_search"/*.a3m; do
    base=$(basename "$file" .a3m)
    echo "Processing $base"

    # Extract the second line (sequence) and remove newline
    seq_line=$(sed -n '2p' "$file" | tr -d '\r\n')
    seq_length=${#seq_line}

    # Generate subunit files
    for s in $(seq 1 "$max_subunits"); do
        out_file="$output_dir/${base}.${s}.a3m"
        {
            echo "# ${seq_length}	${s}"
            echo ">${base}"
            echo "${seq_line}"
            # Also include a cleaned version of the original
            sed -e 's/[[:space:]]\+/ /g' "$file"
        } > "$out_file"
    done
done
