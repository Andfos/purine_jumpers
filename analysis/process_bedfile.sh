#!/bin/bash

# Default values
sorted=false
complement=false
extend=0

# Function to display usage information
usage() {
	echo "Usage: $0 -i genome_index -b bed_file1 -o output_file [-s] [-c] [-e]"
	exit 1
}

# Process command line options
while getopts ":i:b:o:e:s:c" opt; do
	case ${opt} in
		i)
		genome_index=$OPTARG
		echo "Genome index is ${genome_index}."
		;;
		b)
		bed_file=($OPTARG) # use the split+glob operator
		echo "Bedfile is ${bed_file}."
		;;
		o)
		output_file=$OPTARG
		echo "Output file is ${output_file}."
		;;
		e)
		extend=$OPTARG
		echo "Extension is ${extend}."
		;;
		s)
		sorted=true
		;;
		c)
		complement=true
		;;
		\?)
		echo "Invalid option: -$OPTARG" >&2
		usage
		;;
		:)
		echo "Option -$OPTARG requires an argument." >&2
		usage
		;;
	esac
done


# Check for mandatory options
if [ -z "$genome_index" ] || [ -z "$bed_file" ] || [ -z "$output_file" ]; then
  echo "Missing required options."
  usage
fi


# Create a temporary directory and instantiate temporary files.
temp_dir=$(mktemp -d)
temp_sorted_bed="${temp_dir}/temp_sorted.bed" 
temp_sorted_genome_index="${temp_dir}/sorted_genome_index.fai"
temp_slop_bed="${temp_dir}/temp_slop.bed" 


# Sort the files if necessary (sorted=false).
if [ "$sorted" = false ]; then
	gsort -k1,1V -k2,2n $bed_file > $temp_sorted_bed
	gsort -k1,1V -k2,2n $genome_index > $temp_sorted_genome_index
else
	$temp_sorted_bed = $temp_bed
	$temp_sorted_genome_index = $genome_index
fi



# Add +/- <extend> bases to the beginning and end of each range in the bedfile
# (extend != 0).
if [ "$extend" != 0 ]; then
	bedtools slop -i $temp_sorted_bed -g $temp_sorted_genome_index -b $extend > $temp_slop_bed
else
	echo "extend set to 0"
	$temp_slop_bed = $temp_sorted_bed
fi




# Find the complement of the file if necessary (complement=false)
if [ "$complement" = true ]; then
	bedtools complement -i "$temp_slop_bed" -g "$temp_sorted_genome_index" > "$output_file"
else
	cat "$temp_sorted_bed" > "$output_file"
fi


# Echo that the script has finished running.
echo "Bedfile processed. Results saved in $output_file"



# Remove the temporary directory and its contents
rm -r "$temp_dir"


