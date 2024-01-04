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
the intergenic, nucleosome-bound, and highly-conserved region of the *hg19* human
genome for downstream processing.

1. Obtain a genome file and index file. For example, to obtain the *hg19* genome, 
   simply proceed to the following link:
    ```
    https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    ```
   Following download, you can unzip the compressed file using 
   `gunzip <path/to/genome_file>`.

   
   We must next obtain a genome index file. This can be accomplished
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

   a) **Obtain the genomic regions overlapping genes and save to a bed-file**. 
   We proceed first to the UCSC Table Browser located at 
   `https://genome.ucsc.edu/cgi-bin/hgTables`. Here, we specify the 
   following parameters:
            
        clade    : Mammal
        genome   : Human
        assembly : Feb. 2009 (GRCh37/hg19)
        group    : Genes and Gene Predictions
        track    : UCSC Genes
        table    : knownGene
        region   : genome
        output format : BED - browser extensible data
       
   We can then specify the output filename and download the BED file using 
   the `get output` button.

   This will only retrieve the regions containing genes within the *hg19* 
   genome. In order to obtain intergenic regions, we use the
   `process_bedfile` script in the `analysis` directory: 
        
        bash analysis/process_bedfile.sh -i data/hg19.fa.fai -b data/hg19_genes.bed -e 500 -c -o data/hg19_intergenic.bed
        
   In the example, we specify that we will extend each region in the input BED file 
   by 500 nucleotides using the -e option. We provide the -c option, which 
   will obtain the complement of the extended regions of the input BED file. 
   Our output file will therefore contain the intergenic regions of the
   *hg19* genome that are separated from known genes by 500 
   nucleotides.


   b) **Obtain the conserved regions of the genome**. To accomplish this, we will 
   use results from the *PhastCons-100Way* experiment (details can be found 
[here](https://genome.ucsc.edu/cgi-bin/hgc?hgsid=916826631_g8XasCQqrg8t9dxczEQmzhNA9Nyc&c=chr12&l=53858048&r=53859044&o=53858048&t=53859044&g=phastCons100way&i=phastCons100way)). To access this data, we can connect to the table using MySQL via the 
   MariaDB. The MySQL workbench is a good tool for this task, although one 
   could also use the CLI (see [here](http://genome.ucsc.edu/goldenPath/help/mysql.html)
   for details on how to connect to the MariaDB server). 
   Once connected to the server, we can run the following query to retrieve 
   all genomic bins where the sumData column is >= 165: 
            
        SELECT chrom, chromStart, chromEnd
        FROM hg19.phastCons100way
        WHERE sumData >= 165

   There are 2,813,440 rows in the table. By filtering for rows where the sumData 
   is >165, we are left with 560,596 rows. Therefore, we are obtaining approximately 
   the top-20% most conserved 1024-nucleotide stretches of the *hg19* genome. 
   Once the results of the query are returned, export the data as a TSV file 
   in BED format.


   c) **Obtain the nucleosome-bound regions of the genome**. To find the 
   positions of the genome bound by nucleosomes, we will use the results 
   of an experiment from an Mnase digestion of GM12878 cells. The details of 
   the experiment can be found [here](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeSydhNsome&hgta_table=wgEncodeSydhNsomeGm12878Sig&hgta_doSchema=describe+table+schema). Access and download the 
   file titled *wgEncodeSydhNsomeGm12878Sig.bigWig* from the data archive 
   [here](https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhNsome/). 
   
   Once downloaded, we will need to convert the file from bigWig format to 
   BED format. This will entail using another BED file with predefined 
   regions, and then calculating an average coverage score of those regions 
   using the bigWig file. Therefore, we will first create a BED file for the 
   *hg19* genome with non-overlapping 500 nucleotide regions using the 
   `make_bed_from_index.py` script in the `analysis` folder:

        python analysis/make_bed_from_index.py -i data/hg19.fa.fai -o data/hg19_500nt.bed -n 500

   We can then calculate the average nucleosome-signal over these regions using 
   the `bigWigAverageOverBed` UCSC tool, available for download 
   [here](http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/). Once the 
   tool is downloaded and made executable using `chmod +x bigWigAverageOverBed`, 
   use the following command:

        bigWigAverageOverBed data/wgEncodeSydhNsomeGm12878Sig.bigWig data/hg19_500nt.bed -bedOut=data/hg19_nucCoverage.bed temp.tab > /dev/null

   **Note:** To use the UCSC tool without specifying the full-path, you will need 
   to update the shell's PATH variable accordingly. 

   Once we have obtained the `data/hg19_nucCoverage.bed` file with the average 
   nucleosome coverage across all non-overlapping 500 nucleotide regions in 
   *hg19*, we will need to filter it to obtain the regions of high nucleosome 
   occupancy:

        awk -F'\t' '$5 > 1.32' data/hg19_nucCoverage.bed > data/hg19_nucleosome.bed

   The original file contains 6,274,368 regions, and the resulting file
   contains 1,255,658 regions. We have therfore retrieved approximately the 
   top-20% most nucleosome bound regions in the *hg19* genome. 

 

   d) **Find overlap among bed files**. We have now obtained 3 bed files that 
   we must combine to find the regions of the hg19 genome that are intergenic, 
   conserved, and bound by nucleosomes. We can accomplish this using 
   `combine_bedfiles.sh` in the `analysis` directory:

        bash analysis/combine_bedfiles.sh -b data/hg19_intergenic.bed,data/hg19_nucleosome.bed,data/hg19_conserved.bed -o data/hg19_inNucCon.bed

   We obtain a total of 41,113,014 nucleotides in our final sequence,
   corresponding to about 1.4% of the entire *hg19* genome.
 

   e) **Extract a FASTA file from the combined-bedfile.** We must now use our 
   genome file and our final BED file to obtain the actual 
   DNA sequence of our retrieved regions. We do this using the following command:

        bedtools getfasta -fi data/hg19.fa -bed data/hg19_inNucCon.bed -fo data/hg19_inNucCon.fa
   
   f) **Process the sequence.** Once we have our final FASTA file, we can then 
   process it as illustrated above in the ***Running the Analysis*** section.
                  







