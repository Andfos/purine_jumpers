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
echo "Reults written to $outfile. All sequence data is on one line."

# Recode purins to R and pyrimidines to Y. Write all results to a single line of
# the output file.
sed '/^[^>]/s/[ag]/R/gI' $infile | sed '/^[^>]/s/[tc]/Y/gI' | 
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > $outfile
