#!/bin/bash


# Parse command line arguments.
#while getopts ":i:o:" opt; do
#  case $opt in
#    i)
#      infile="$OPTARG"
#      ;;
#    o)
#      outfile="$OPTARG"
#      ;;
#    \?)
#      echo "Invalid option: -$OPTARG" >&2
#      exit 1
#      ;;
#    :)
#      echo "Option -$OPTARG requires an argument." >&2
#      exit 1
#      ;;
#  esac
#done
#
#
# Explain to user the progress.
#echo "Recoding $infile"
#echo "Results written to $outfile. All sequence data is on one line."

# Recode purins to R and pyrimidines to Y. Write all results to a single line of
# the output file.
#sed '/^[^>]/s/[ag]/R/gI' $infile | sed '/^[^>]/s/[tc]/Y/gI' | 
#	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > $outfile

#!/bin/bash

# Parse command line arguments.
while getopts ":i:o:" opt; do
  case $opt in
    i)
      infile="$OPTARG"
      ;;
    o)
      outfile="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Explain to user the progress.
echo "Recoding $infile"
echo "Results written to $outfile. Each sequence data is on its own line."

# Process each region separately, recode, and write to the output file.
awk '/^>/ { if (seq) print seq; printf("%s\t", $0); seq=""; next } { seq = seq $0 } END { print seq }' "$infile" |
while read -r region; do
    # Extract header and sequence from the region
    header=$(echo "$region" | awk -F'\t' '{print $1}')
    sequence=$(echo "$region" | awk -F'\t' '{print $2}')

    # Recode purines to R and pyrimidines to Y
    recoded_sequence=$(echo "$sequence" | sed 's/[ag]/R/gI; s/[tc]/Y/gI')

    # Output the recoded region to the output file
    echo -e "$header\n$recoded_sequence" >> "$outfile"
done

echo "Processing complete. Output written to '$outfile'."

