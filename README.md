# kmer_id
metagenomic read identification by kmer database

source file newkmer_10nx.cpp
uses gzip library, compile with:
g++ -O3 newkmer_10nx.cpp -o nk10 -lz

database files needed in bact10 subdirectory:
bData10.txt
btree10.txt
refkey10.txt
probes10.txt.gz - too large (1.5 Gb) for github, email me for box link.

run (uses 25 Gb RAM) with
./nk10 /path-to-fastq-files/
Input files are two paired trimmed fastq files (_R1_tr.fastq.gz and _R2_tr.fastq.gz - names can be changed in source)

Python 3 program readbatch_10.py  
to collect results files and generate report
Change dir1 to location of result files
Output is .csv file with read count and % abundance for each species.
