# WPS
Calculate average WPS for genes

Imput files:
1. WPS.bw
2. human_gene.txt

Command lines:
convert bigwig to bedgraph: -bigWigToBedGraph wps.bw wps.bedGraph

Python algorithms:
1. pandas read files (wps.bedGraph, human_gene.txt)
2. groupby as chr
3. loops the chromosomes
4. loops the genes in each chromosome
5. loops the wps:
    a. change the starting/ending position to the positions of genes
    b. the WPS with starting position <= ending position are those we are looking for inside the gene
    c. calculate the avg. WPS and store in the dictionary(key: gene version ID, value: avg. WPS)

There are two .py files under this repository, bioinfo_test.py is for testing (without looping chr part) and bioinfo.py is the complete version.

comments:
the whole WPS file is too big for uploading
