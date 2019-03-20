# kmer_id
Mitichondrial read identification by kmer database for Galaxy server

File kmer_read_m3.cpp
uses gzip library, compile with:
g++ -O3 -std=c++0x kmer_read_m3.cpp -o kmerread -lz

database files needed in same directory:
(can create with kmer_build, but not described here yet)
mitochondria_data.txt
mitochondria_refkey.txt
mitochondria_count.txt
mitochondria_tree.txt
mitochondria_probes.txt.gz
1a.fasta (test input)

File kmer_read_m3.py  
(Python 2.7)
Run with:
python kmer_read_m3.py -w [working directory] -d [output directory] -i [input filename1] [input filename2]

Input files can be two paired files (.fastq, .fastq.gz, .fasta, .fasta.gz) or a single file with none as filename2
Output is .csv file with read count and % abundance for each species.
