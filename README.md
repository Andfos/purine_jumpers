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
2. Next, invoke `python analysis/find_jumpers.py -i <path_to_input_file>
   -o <path_to_output_file>`, along with any other command line arguments that
the user wishes to change from the default (run `python analysis/find_jumpers.py
-h` to display a list of all command line arguments). This script will move
along the purine-coded fasta file one nucleotide at a time, and search for
regions of homology between the first `<reference_size>; default = 40` nucleotides and the next
`<window_size - reference_size>` nucleotides. If a single region of `<reference_size>` 
is greater than or equal to `<min_homology>; default = 0.6` and has a unique, highest homology within
the window, it will be deemed a jumper and written to the output file.
3. Next, invoke `find_odd_even_homology.py -i <path_to_input_file> -o
   <path_to_output_file>`. This script will split the jumpers into phases and
calculate the difference in homology between the odd-numbered phase and
the even-numbered phase. 
4. Lastly, the `jumpers_analysis.ipynb` can be run to generate plots to
   visualize the results of the experiment.

## Example of use
In the following section, we provide an example of how to run the code on a
randomly-generated FASTA file consisting of 1000 nucleotides. Note that at the
second position of the file, there is an N. This is a masked base, and the
`find_jumpers.py` script will ignore any windows containing an N. 

1. Preprocess the file:
    ```
    bash analysis/preprocessing.sh -i test/dna_1000nt.fa  -o test/dna_1000nt_oneLine.fa
    ```

2. Find jumpers:
    ```
    python analysis/find_jumpers.py -i test/dna_1000nt_oneLine.fa -o test/dna_1000nt_homology_results.tsv -w 200 -r 40 -m 0.60
    ```

3. Find odd-even phase homology:
    ```
    python analysis/find_odd_even_homology.py -i test/dna_1000nt_homology_results.tsv -o test/dna_1000nt_homology_results_oddEven.tsv
    ```

4. Invoke jupyter notebook, and run the code cells in `jumpers_analysis.ipynb`.
   Note that the results displayed on Github are currently those for Human chromosome 20 of HG38 genome. 


## Obtaining custom genomic regions
In this section, we explain how one can obtain customized genomic regions for 
analysis with this package. We illustrate this with an example where we obtain 
the intergenic, nucleosome-bound, and highly-conserved region of the humand
genome for downstream processing.

1. Obtain a genome file and index file. For example, to obtain the *hg19* genome, 
   simply proceed to the following link:
    ```
    https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    ```
   Following download, obtain a genome index file. This can be accomplished
   using either of the following two commands:
    ```
    bioawk -c fastx '{print $name"\t"length($seq)}' <genome_file> | gsort -k1,1V -k2,2n)
    ```
   or
    ```
    samtools faidx <genome_file>
    ```

2. Obtain required bedfiles. Useful tools include the UCSC Table Browser and 
   the MariaDB.

    a) Obtain the intergenic regions. 
       We proceed first to the UCSC Table Browser located at
       `https://genome.ucsc.edu/cgi-bin/hgTables`. Here, we specify the 
       following parameters:
        ```
        clade    : Mammal
        genome   : Human
        assembly : Feb. 2009 (GRCh37/hg19)
        group    : Genes and Gene Predictions
        track    : UCSC Genes
        table    : knownGene
        region   : genome
        output format : BED - browser extensible data
        ```
       We can then specify the output filename and download the BED file using 
       the `get output` button.

