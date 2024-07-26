#!/bin/bash

max_subunits=10

# Input directory file
input_search=$1

# Output directory
output_dir=$2

# create output directory 
mkdir -p $output_dir

# for each file in the input directory
for i in `ls $input_search`; do
	echo $i
    for s in $(seq 1 $max_subunits); do

        # extract the length of the sequence
	seq_length=$(sed -n '2p' "$input_search/$i/0.a3m" | wc -c)
        seq_length=$((seq_length - 1))

        # create the output file with the desired content
        { echo "# $seq_length	$s"; echo "> $i"; sed -n '2p' "$input_search/$i/0.a3m" ; sed -e 's/[[:space:]]\+/ /g'  "$input_search/$i/0.a3m"; } > "$output_dir/$(basename "${i%.*}").$s.a3m"
    done
done

