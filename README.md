This code allows a user to search for regions of homology nearby in a
chromosome, as well as to analyze the differences in homology between even and
odd phases of nucleotides.

## Running the analyis
1. A fasta file containing a single DNA sequence can be purine coded and placed
 on a single line by invoking `analysis/preprocessing -i <path_to_input_file> -o
<path_to_output_file>`. This will recode all adenine (A) and guanine (G)
nucleotides to R, and all thymine (T) and cytosine (C) to Y. It will also place
the entire sequence on a single line, which is useful for the next step in the
analysis.
2. Next, the user should invoke `python analysis/find_jumpers.py -i <path_to_input_file>
   -o <path_to_output_file>`, along with any other command line arguments that
the user wishes to change from the default.
