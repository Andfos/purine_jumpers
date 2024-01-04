#!/bin/bash

# Default values
sorted=false

# Function to display usage information
usage() {
	echo "Usage: $0 -b bed_file1 [bed_file2 bed_file3 ...] -o output_file [-s]"
	exit 1
}

# Process command line options
while getopts ":b:o:s" opt; do
	case ${opt} in
		b)
		set -f # disable glob
		IFS=',' # split on comma characters
		bed_files=($OPTARG) # use the split+glob operator
		;;
		o)
		output_file=$OPTARG
		;;
		s)
		sorted=true
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
if [ -z "$bed_files" ] || [ -z "$output_file" ]; then
	echo "Missing required options."
	usage
fi


# Make a temporary directory to store intermediate files.
temp_dir=$(mktemp -d)


# Sort input bed files if not already sorted
if [ "$sorted" = false ]; then
    for bed_file in "${bed_files[@]}"; do
        #sorted_bed_file="${temp_dir}/${bed_file%.bed}_sorted.bed"
		
		# Use basename to extract the filename
		filename=$(basename "$bed_file")
		sorted_bed_file="${temp_dir}/${filename%.bed}_sorted.bed"
        
        # Ensure the directory for the sorted files exists
        mkdir -p "$temp_dir"
        
        # Use gsort to sort the bed file and save it to the temporary directory
        #gsort -k1,1V -k2,2n "$bed_file" > "${sorted_bed_file}"
        gsort -k1,1V -k2,2n "$bed_file" | awk -F'\t' '$1 ~ /^(chr[1-9]|chr1[0-9]|chr[XY])$/' > "${sorted_bed_file}"

        # Store the path to the sorted bed file
        sorted_bed_files+=("${sorted_bed_file}")
    done
else
    sorted_bed_files=("${bed_files[@]}")
fi

# Print the paths of the sorted bed files
#printf "%s\n" "${sorted_bed_files[@]}"




# Combine the bedfiles.
combined_bed="${temp_dir}/combined.bed"
#multiIntersectBed -i "${sorted_bed_files[@]}" > "${combined_bed}"
bedtools multiinter -i "${sorted_bed_files[@]}" > "${combined_bed}"


# Extract rows where all bedfiles overlap. 
num_columns=$(awk '{print NF; exit}' "$combined_bed")
num_bedfiles=$((num_columns - 5))
awk -v num_bedfiles="$num_bedfiles" -F'\t' '$4 == num_bedfiles' "${combined_bed}" > "${output_file}"



# Use awk to calculate the sum of the start positions and end positions.
sum_start=$(awk '{ sum += $2 } END { print sum }' "$output_file")
sum_end=$(awk '{ sum += $3 } END { print sum }' "$output_file")


# Calculate the total sequence length of the final bedfile.
total_sequence_length=$((sum_end - sum_start))
echo "Total sequence length: ${total_sequence_length}"





# Remove the temporary directory and its contents
rm -r "$temp_dir"



