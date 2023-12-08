import argparse


# Parse CLI arguments.
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True, help="Input file path")
    parser.add_argument("-o", "--output_file", required=True, help='Output file path')
    parser.add_argument(
            "-w", "--window_size", default=200, 
            help="Window size for nucelotides. Homology will be searched in the region <window_size> - <reference_size> (default: 200)")
    parser.add_argument(
            "-r", "--reference_size", default=40, 
            help="Number of nucleotides to to use as reference for finding jumpers (default: 40)")
    parser.add_argument(
            "-m", "--min_homology", default=0.6, 
            help="Minimum homology value to be considered a jumper (default : 0.6)")

    args = parser.parse_args()
    return args



# Find the best matching string to a reference sequence.
def find_best_matching_string(reference, sequence, min_homology=0.6):
    
    # Initialize variables.
    best_homology = 0
    best_match_string = ""
    julen = 0
    homology_dict = {}
    
    # Iterate along sequence.
    for i in range(len(sequence) - len(reference) + 1):
        substring = sequence[i:i + len(reference)]

        # Count the number of matching characters and obtain homology.
        match_count = sum(c1 == c2 for c1, c2 in zip(reference, substring))
        homology = match_count / len(reference)
        
        # If the homology is greater than the best homology, record info.
        if (homology >= best_homology):
            best_homology = homology
            jumper = substring
            julen = i + len(reference)
            
            # Add information to a homology dictionary with keys as homology
            # values.
            if homology not in homology_dict:
                homology_dict[homology] = []
            temp_dict = {
                    "homology" : homology,
                    "jumper" : jumper,
                    "julen" : julen}
            homology_dict[homology].append(temp_dict)
                                   
    return homology_dict





# Process the input file.
def process_file(
        input_file_path, output_file_path, 
        window_size, reference_size, min_homology):
    

    with open(input_file_path, "r") as file, open(output_file_path, "w") as output_file:
        
        # Write the header of the output file.
        output_file.write("position\thomology\tjulen\treference\tjumper\n")
        
        # Read the first 200 characters
        window = file.read(window_size)
        
        position = 0
        while len(window) == window_size:
            nt = window[0]
            
            # Skip past windows that contain N
            if "n" in window.lower():
                window = window[1:] + file.read(1)
                position += 1
                continue
            
            # Monitor progress by printing position every 100,000 nt.
            if (position % 100000) == 0:
                print(position)

            # Take the first 40 characters as the reference
            reference = window[:reference_size]

            # Find the best matching string within the next 160 characters
            homology_dict = find_best_matching_string(
                    reference, window[reference_size:], min_homology)
            
            
            # Only retrieve jumpers that have max homology, and only if there 
            # are no two or more jumpers with the same max homology.
            max_homo = max(homology_dict.keys())
            if (len(homology_dict[max_homo]) > 1) or (max_homo <= min_homology):
                position += 1
                window = window[1:] + file.read(1)
                continue
            
            # Retrieve homology, jump length (julen), and jumper string.
            homology = homology_dict[max_homo][0]["homology"]
            julen = homology_dict[max_homo][0]["julen"]
            jumper = homology_dict[max_homo][0]["jumper"]

            # Write results to file
            output_file.write(
                    f"{position}\t{homology}\t{julen}\t{reference}\t{jumper}\n"
            )
            
            # Slide the window one character to the right
            window = window[1:] + file.read(1)
            position += 1





# Main function to find jumpers.
def main():

    # Parse CLI arguments.
    args = parse_arguments()
    input_file_path = args.input_file
    output_file_path = args.output_file
    window_size = int(args.window_size)
    reference_size = int(args.reference_size)
    min_homology = float(args.min_homology)
    
    # Process the file to find jumpers.
    process_file(
            input_file_path, output_file_path, window_size, reference_size,
            min_homology)




if __name__ == "__main__":
    main()
