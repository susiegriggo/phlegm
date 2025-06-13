# Choosing the best representative - one with the highest bitscore to all

```bash
FILE="per_phrog_tsvs/phrog_2.tsv.gz"
gunzip -c $FILE | awk -F'\t' '{sum[$1] += $12} END {for (key in sum) print key, sum[key]}' | sort -k2,2n | tail -n 1
```

```bash
#!/bin/bash

OUTPUT_FILE="susie_representative_summary.tsv"
echo -e "File\tMax_Key\tMax_Sum" > $OUTPUT_FILE

for i in $(seq 1 38880); do
    FILE="per_phrog_tsvs/phrog_${i}.tsv.gz"
    if [[ -f $FILE ]]; then
        result=$(gunzip -c "$FILE" | grep -v -e "envhog" -e "efam" -e "dgr" | \
                 awk -F'\t' '{sum[$1] += $12} END {for (key in sum) print key, sum[key]}' | \
                 sort -k2,2n | tail -n 1)
        
        if [[ -n "$result" ]]; then
            echo -e "${FILE}\t${result}" >> $OUTPUT_FILE
        fi
    else
        echo -e "${FILE}\tNA\tNA" >> $OUTPUT_FILE
    fi
done

echo "Processing complete. Results saved in $OUTPUT_FILE."

```

```bash
OUTPUT_FILE="pharokka_new_representative_summary.tsv"
echo -e "File\tMax_Key\tMax_Sum" > $OUTPUT_FILE

for i in $(seq 1 38880); do
    FILE="per_phrog_tsvs/phrog_${i}.tsv.gz"
    if [[ -f $FILE ]]; then
        result=$(gunzip -c "$FILE" |  \
                 awk -F'\t' '{sum[$1] += $12} END {for (key in sum) print key, sum[key]}' | \
                 sort -k2,2n | tail -n 1)
        
        if [[ -n "$result" ]]; then
            echo -e "${FILE}\t${result}" >> $OUTPUT_FILE
        fi
    else
        echo -e "${FILE}\tNA\tNA" >> $OUTPUT_FILE
    fi
done

echo "Processing complete. Results saved in $OUTPUT_FILE."
