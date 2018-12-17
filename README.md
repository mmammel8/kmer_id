# kmer_id
Mitichondrial read identification by kmer database

File kmer_read_m1.cpp
uses gzip library, compile with:
g++ -O3 -std=c++0x kmer_read_m1.cpp -o kmerread -lz

database files needed in subdirectory /mitochondria:
(supplied by me or can create with kmer_build, but not described here yet)
mitochondria_data.txt
mitochondria_key.txt
mitochondria_count.txt
mitochondria_tree.txt
mitochondria_probes.txt

File kmer_readm.py  
(Python 2.7)
Input file: jobs.txt
(actual filename is specified in line 9 to jobs_name, but should be changed to a command line argument)

format of input file:
<p>sample_name1 number_of_input_files1
<p>inputfile_1-1
<p>inputfile_1-2
<p>sample_name2 number_of_input_files2
<p>inputfile_2-1
<p>inputfile_2-2
<p>...

example:
<p>mRNA24hC1 2
<p>/mnt/gno/gnome3/DropBox/mRNA24hC1_S1_L001_R1_001.fastq.gz
<p>/mnt/gno/gnome3/DropBox/mRNA24hC1_S1_L001_R2_001.fastq.gz
<p>mRNA24hC2 2
<p>/mnt/gno/gnome3/DropBox/mRNA24hC2_S2_L001_R1_001.fastq.gz
<p>/mnt/gno/gnome3/DropBox/mRNA24hC2_S2_L001_R2_001.fastq.gz

Output is .csv file with read count and % abundance for each sample.
