import argparse


# Parse CLI arguments.
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--index", required=True, help="Genome index file path")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument(
            "-n", "--step", default=500, 
            help="Step size for creating genomic regions (default: 500)")
    args = parser.parse_args()
    return args



# Open input and output files.
def write_bed(index, output, step):
    with open(index, "r") as index_file, open(output, "w") as output_file:
        
        # Initialize region numbers to 1.
        region_num = 1
        for line_number, line in enumerate(index_file, 1):
            
            # Assuming the first two columns of index file are: chromosome_name, chromosome_length
            fields = line.strip().split("\t")
            chr_name = fields[0]
            chr_length = int(fields[1])

            # Create non-overlapping BED entries for every <step> nucleotides in the genome.
            for current_start in range(0, chr_length, step):
                current_end = min(current_start + step, chr_length)
                bed_entry = (
                        f"{chr_name}\t{current_start}\t{current_end}\t{region_num}\n"
                )
                output_file.write(bed_entry)
                region_num += 1


# Main function.
def main():
    args = parse_arguments()
    index =  args.index
    output = args.output
    step = int(args.step)
    write_bed(index, output, step)



if __name__ == "__main__":
    main()
