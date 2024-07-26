#!/bin/bash

# Input FASTA file
input_fasta=$1

# Output directory
output_dir=$2
mkdir -p $output_dir

# Read the input FASTA file
while read -r line; do
    if [[ $line == ">"* ]]; then
        # Remove the leading '>' and any spaces from the title
        title=$(echo "$line" | sed 's/^>//' | tr -d '[:space:]')
        # Create a new file with the sequence title
        file_name="${output_dir}/${title}.fasta"
        echo "$line" > "$file_name"
    else
        # Append the sequence to the current file
        echo "$line" >> "$file_name"
    fi
done < "$input_fasta"

echo "FASTA entries have been separated into individual files in the '$output_dir' directory."

