// newkmer match k-mers to reads
// opt verify reads by SmithWaterman
// MKM 1 APR 2016
//-std=c++0x -D__NO_INLINE__   ... at end: -lz
//#include "stdafx.h"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <time.h>
//#include <sys/time.h>
#include <cstddef>
#include <cstring>
#include <algorithm>
#include <set>
#include <zlib.h>
#include <dirent.h>
using namespace std;

const int minalign = 0; //number of reads to align per species
#define FASTQ 1  //fasta = 0
string e1 = "_R1_tr.fastq.gz";
string e2 = "_R2_tr.fastq.gz";
//string e1 = "_L001_R1_001.fastq.gz";
//string e2 = "_L001_R2_001.fastq.gz";		
//string e1 = "_R1.fasta";

typedef int ptype;
typedef unsigned long long ui64;
typedef unsigned long long ktype;
typedef unsigned int vtype;
typedef unsigned short int otype;
typedef unsigned long long itype;

//int TARGET = 0; //save these, second input parameter
const int KSIZE = 30; //kmer size
const int MAXORGS = 14791;
const int MAXTAR = 5982; //2 to MAXTAR-1 inclusive are present, 1 is root, 0 unused
const int MAXREP = 2048;
const int REPSHIFT = 11; //(log of maxrep)
const int SAVENUM = 12; //save this many reads
const itype MAXHASH = (1 << 30); //pow2
//const itype MAXHASH = 1048576*64;
const int GAPO =11;
const int GAPX =1;
const int MATCH = 5;
const int MISMATCH = -4;
const bool INIGAPPEN = false;
//const int threshold_scr = 900;
const int MAXALIGN = 1024;
const int beam = 8;
string accession[MAXORGS];
int targno[MAXORGS], ntargorgs[MAXTAR], mcount;
int gcount[MAXTAR], savepos, tct, norgs=0;
int ucount[MAXTAR];
ofstream outread;
set<ktype> kmer_seen;

//input data files
string iname = "./bact10/bData10.txt";   //text file of groups for probe design
string pname = "./bact10/probes10.txt.gz";
string tname = "./bact10/btree_10.txt";   //edges in taxonomy tree
string dir = "/mnt/mgb/Mark_backup/genbank/"; // directory for reference sequences for alignment

const char bases[4] = { 'A', 'C', 'G', 'T' };
const ui64 one = 1;
const ui64 two = 2;
const ui64 three = 3;
const ui64 mask = (one << (KSIZE * 2)) - 1;
const ui64 hiC = one << (KSIZE - 1) * 2;
const ui64 hiG = two << (KSIZE - 1) * 2;
const ui64 hiT = three << (KSIZE - 1) * 2;
const ui64 hi1 = one << 62;
const ui64 hi2 = two << 62;
const ui64 hi3 = three << 62;
const ui64 himask = hi1 - one;

const unsigned BUFLEN = 0x4000;

void error(const char* const msg)
{
    cerr << msg << "\n";
    exit(255);
}

class Tree1
{
public:

	vector<int> children[MAXTAR];
	int parent[MAXTAR];
	int root;

	Tree1()
	{
		root = 1;
		for (int i = 0; i < MAXTAR; i++)
			parent[i] = root;
	}

	~Tree1()
	{
	}

	void add_edge(int x, int y)
	{
		parent[y] = x;
		children[x].push_back(y);
	}

	int msca(int x, int y)
	{
		//return most specific common ancestor of nodes x and y
		set<int> ancestors;
		int z;

		ancestors.insert(root);
		z = x;
		//cout << x << " " << y << " " << z << endl;
		while (z != root)
		{
			ancestors.insert(z);
			z = get_parent(z);
		//cout << x << " " << y << " " << z << endl;
		}
		z = y;
		if (ancestors.find(y) != ancestors.end())
			return x; //x is descendent of y
		while (ancestors.find(z) == ancestors.end())
		{
			z = get_parent(z);
			if (z==x) //y is descendent of x
				return y;
		//cout << x << " " << y << " " << z << endl;
		}
		return z; //common ancesctor
	}

	int get_parent(int x)
	{
		if (x != root && x > 0)
			return parent[x];
		else
			return root;
	}

};

Tree1 *taxonomy;

class Hashtable
{
public:

	itype size;

	struct Cell
	{
		ktype key;
		vtype value;
		int org;
		ptype position;
 		bool fstrand;
	} *pHash;

	Hashtable()
	{
		size = 0;
		pHash = new Cell[MAXHASH];
		HashClear();
		//for (int i = 1; i<MAXREPROBE; i++)
		//	reprobeVal[i] = reprobeVal[i - 1] + i; //i*i/2+i/2
	}

	~Hashtable()
	{
		delete[] pHash;
		size = 0;
	}

	// from code.google.com/p/smhasher/wiki/MurmurHash3
	itype integerHash(ui64 k)
	{
		k ^= k >> 33;
		k *= 0xff51afd7ed558ccd;
		k ^= k >> 33;
		k *= 0xc4ceb9fe1a85ec53;
		k ^= k >> 33;
		return (itype) k;
	}

	void HashClear()
	{
		memset(pHash, 0, MAXHASH * sizeof(Cell));
	}

	int getHash(ktype key, otype &org, ptype &position, bool &fstrand)
	{ 
		//0 target=bad
		//returns target if probe present in table
		//sets position and fstrand
		itype index, reprobe, i, hash;

  		position = 0;
		hash = integerHash(key);
		reprobe = 0;
		i = 0;
		do
		{
			index = ((hash + reprobe) & (MAXHASH-1)); //for power of 2 size
			reprobe += ++i;
			if (pHash[index].value == 0)
			{ //empty cell
				return 0;
			}
			else if (pHash[index].key == key)
			{
				position = pHash[index].position;
				fstrand = pHash[index].fstrand;
				org = pHash[index].org;
				return (pHash[index].value);
			}
		} while (reprobe < MAXHASH);

		return 0;
	}

	void add_kmer(ktype key, vtype target, otype org, ptype position, bool fstrand)
	{ //add kmer to table
		itype index, reprobe, i, hash;
		bool notdone=true;

		hash = integerHash(key);
		reprobe = 0;
		i = 0;
		do
		{
			index = ((hash + reprobe) & (MAXHASH-1)); //for power of 2 size
			reprobe += ++i;

			if (pHash[index].value == 0)
			{ //empty cell
				notdone = false;
				pHash[index].key = key;
				pHash[index].value = target;
				pHash[index].org = org;
				pHash[index].position = position;
				pHash[index].fstrand = fstrand;
				if (++size > MAXHASH - 32)
				{
					cout << "out of memory in table " << endl;
					exit(1);
				}
			}
		} while (notdone);
	}

}; //class Hashtable
Hashtable *ht;

string process_seq(string seq)
{  
	//format nucleotide sequence to all caps ACTG, all else is N
	string seq2;
	string::iterator it1;
	char base;

	seq2 = "";
	for (it1 = seq.begin(); it1 != seq.end(); ++it1)
	{
		base = *it1;
		switch (base)
		{
		case 'a':
			seq2 += "A";
			break;
		case 'c':
			seq2 += "C";
			break;
		case 'g':
			seq2 += "G";
			break;
		case 't':
			seq2 += "T";
			break;
		case 'A':
			seq2 += "A";
			break;
		case 'C':
			seq2 += "C";
			break;
		case 'G':
			seq2 += "G";
			break;
		case 'T':
			seq2 += "T";
			break;
		default:
			seq2 += "N";
			break;
		}
	}

	return seq2;
}

string load_data3(gzFile in)
{ //read gzip fasta file and extract nt sequence
	char buf[BUFLEN];
	char* offset = buf;
	string line;
        string sequence = "";

	while (true) 
	{
		int err, len = sizeof(buf)-(offset-buf);
		if (len == 0) error("Buffer to small for input line lengths");

 		len = gzread(in, offset, len);

		if (len == 0) break;    
		if (len <  0) error(gzerror(in, &err));

		char* cur = buf;
		char* end = offset+len;

		for (char* eol; (cur<end) && (eol = find(cur, end, '\n')) < end; cur = eol + 1)
 		{
			line = string(cur, eol);
			if (line.back() == '\r')
				line.pop_back();
			if (line.length() > 0)
			{
				if (line.at(0) == '>')
				{
					sequence += "N"; //contig separator
				}
				else
				{
					sequence += process_seq(line);
				}				
			}
			//cout << line << endl;
		}

		// any trailing data in [eol, end) now is a partial line
		offset = copy(cur, end, buf);
	}

	// trailing data without eol
	line = string(buf, offset);
	if (gzclose(in) != Z_OK) error("failed gzclose");

	return sequence;
}


int *tableM, *tableI, *tableD;

void tableinit()
{
	int i,j;
	const int NINF = -2000000000;

	tableM[0] = 0;
	tableI[0] = NINF;
	tableD[0] = NINF;
	for (i = 1; i<MAXALIGN; i++) //initialize upper row
	{
		j = (0 > i - beam - 1) ? 0 : i - beam - 1;
		tableI[j*MAXALIGN  + i] = NINF;
		if (INIGAPPEN)
		{
			tableD[j*MAXALIGN  + i] = -GAPO - (i - 1)*GAPX;
			tableM[j*MAXALIGN  + i] = -GAPO - (i - 1)*GAPX;
		}
		else
		{
			tableD[0 + i] = 0;
			tableM[0 + i] = 0;
		}
	}
	for (j = 1; j<MAXALIGN; j++) //init left column
	{
		i = (0 > j - beam - 1) ? 0 : j - beam - 1;
		tableD[j*MAXALIGN + i] = NINF;
		if (INIGAPPEN)
		{
			tableI[j*MAXALIGN + i] = -GAPO - (j - 1)*GAPX;
			tableM[j*MAXALIGN + i] = -GAPO - (j - 1)*GAPX;
		}
		else
		{
			tableI[j*MAXALIGN + i] = 0;
			tableM[j*MAXALIGN + i] = 0;
		}
	}
}

int align(string dna1, string dna2)
{
	//smith waterman alignment with beam
	//returns alignment score
	int b1, b2, len1, len2, tindex;
	int b,bb,c,cc;
	int i, j, maxval, startx, starty;
	int i1, i2;

	len1 = dna1.length() + 1;
	len2 = dna2.length() + 1;

	j=1;
	do
	{
		i = (j <= beam) ? j : j - beam;
		i2 = (j + beam > len1) ? len1 : j + beam;
		do
		{
			maxval = max(max(tableM[(j-1)*MAXALIGN + i-1], tableI[(j-1)*MAXALIGN + i-1]), tableD[(j-1)*MAXALIGN + i-1]);
			//copy/replace
                        if (dna1[i-1] != dna2[j-1])
                            maxval += MISMATCH;
			else maxval += MATCH;
			tableM[j*MAXALIGN + i] = maxval;

			b=tableM[(j-1)*MAXALIGN + i] - GAPO; //new gap
			bb = tableI[(j-1)*MAXALIGN + i] -GAPX; //extend
			tableI[j*MAXALIGN + i] = (bb<b) ? b : bb;

			c=tableM[j*MAXALIGN + i-1] - GAPO; //new gap
			cc = tableD[j*MAXALIGN + i-1] -GAPX; //extend
			tableD[j*MAXALIGN + i] = (cc<c) ? c : cc;	
		} while(++i < i2);
	} while(++j < len2);

	i = len1 - 1;
	j = len2 - 1;   //process last entry

	maxval = max(max(tableM[j*MAXALIGN + i], tableI[j*MAXALIGN + i]), tableD[j*MAXALIGN + i]);

	return maxval;

}

int process_read(string &sequence, string acc, int start, int stop)
{
	//returns id of read, 0 = no match, updates gcount, ucount
	//will do align for first minalign of each target
	//will save first SAVENUM reads of each target
	char base;
	int it1, it2;
	string kmer, sequence2, dna1, dna1r, dna2, fname;
	ktype keyF, keyR, key;
	int i, cpos, minscr, readlen2;
	bool fstrand1, fstrand2;
	int readlength = stop - start + 1;
	int st2, scr, stlen2;
	ktype target, final_targ=0;
	ptype position;
	otype org;
	int maxct, totct, rtarg=0;
	bool reject = false;

	cpos = 0;
	keyF = 0;
	keyR = 0;
	minscr = 5 * sequence.length() / 2;
	for (it1 = start; it1 <= stop; ++it1)
	{
		base = sequence.at(it1);
		switch (base)
		{
		case 'A':
			keyF = (keyF << 2) & mask;
			keyR = (keyR >> 2) | hiT;
			cpos++;
			break;
		case 'C':
			keyF = ((keyF << 2) & mask) | 1;
			keyR = (keyR >> 2) | hiG;
			cpos++;
			break;
		case 'G':
			keyF = ((keyF << 2) & mask) | 2;
			keyR = (keyR >> 2) | hiC;
			cpos++;
			break;
		case 'T':
			keyF = ((keyF << 2) & mask) | 3;
			keyR = keyR >> 2;
			cpos++;
			break;
		case 'a':
			keyF = (keyF << 2) & mask;
			keyR = (keyR >> 2) | hiT;
			cpos++;
			break;
		case 'c':
			keyF = ((keyF << 2) & mask) | 1;
			keyR = (keyR >> 2) | hiG;
			cpos++;
			break;
		case 'g':
			keyF = ((keyF << 2) & mask) | 2;
			keyR = (keyR >> 2) | hiC;
			cpos++;
			break;
		case 't':
			keyF = ((keyF << 2) & mask) | 3;
			keyR = keyR >> 2;
			cpos++;
			break;
		default:
			cpos = 0; //ambiguous character
			keyF = 0;
			keyR = 0;
			break;
		}
		if (cpos == KSIZE)
		{
			key = keyF < keyR ? keyF : keyR;
			target = ht->getHash(key, org, position, fstrand2);
			if (target > 0 && gcount[target] < minalign && target != final_targ)
			{ //do alignment on initial reads of a species
				fname = dir + accession[org] + ".fasta.gz";
				sequence2 = load_data3(gzopen(fname.c_str(), "rb"));
				stlen2 = sequence2.length();
				fstrand1 = (keyF < keyR);
				readlen2 = readlength;
				if (fstrand1 == fstrand2)
				{
					st2 = position - it1;
					if (st2 < 0) st2 = 0;
					if (st2 + readlen2 > stlen2) readlen2 = stlen2 - st2;
					dna1 = sequence.substr(start, readlength);
					dna2 = sequence2.substr(st2, readlen2);
				}
				else
				{
					//reverse complement dna1
					st2 = position - KSIZE + 2 + it1 - readlength;
					if (st2 < 0) st2 = 0;
					if (st2 + readlen2 > stlen2) readlen2 = stlen2 - st2;
					//dna1 = sequence.substr(start, readlength);
					dna1 = "";
					for (it2 = stop; it2 >= start; --it2)
					{
						base = sequence.at(it2);
						switch (base)
						{
						case 'A':
							dna1.push_back('T');
							break;
						case 'C':
							dna1.push_back('G');
							break;
						case 'G':
							dna1.push_back('C');
							break;
						case 'T':
							dna1.push_back('A');
							break;
						default:
							dna1.push_back('N');
							break;
						}
					}
					dna2 = sequence2.substr(st2, readlen2);
				}
				scr = align(dna1, dna2);
				//cout << scr << "/" << minscr << endl;
				if (scr < minscr)
				{
					rtarg = target;
					target = 0; //fail
					reject = true;
					//cout << rtarg << ": failed " << scr << endl;
				}

			}
			if (final_targ > 0 && target > 0)
			{
				final_targ = taxonomy->msca(target, final_targ);
			}
			else if (target > 0)
			{
				final_targ = target;
			}
			if (target > 1)
			{
				if (kmer_seen.find(key) == kmer_seen.end())
				{
					ucount[target]++;
					kmer_seen.insert(key);
				}
			}
			cpos--;
		}
	} //it

        if (final_targ > 1 && gcount[final_targ] < SAVENUM)
	{
		//string s = to_string(savepos);
 		outread << ">" << final_targ << ":" << acc << endl << sequence.substr(start, stop-start+1) << endl;
	}
	gcount[final_targ]++;
	tct++; 

	return final_targ;
}

void process_kmer(string sequence, vtype target, otype org, ptype position, bool fstrand)
{	//adds kmer from datafile into hashtable
	char base;
	string::iterator it1;
	int cpos;
	ui64 keyF;

	cpos = 0;
	keyF = 0;
	for (it1 = sequence.begin(); it1 != sequence.end(); ++it1)
	{
		base = *it1;
		switch (base)
		{
		case 'A':
			keyF = (keyF << 2) & mask;
			cpos++;
			break;
		case 'C':
			keyF = ((keyF << 2) & mask) | 1;
			cpos++;
			break;
		case 'G':
			keyF = ((keyF << 2) & mask) | 2;
			cpos++;
			break;
		case 'T':
			keyF = ((keyF << 2) & mask) | 3;
			cpos++;
			break;
		default:
			cpos = 0; //ambiguous character
			keyF = 0;
			break;
		}
		if (cpos == KSIZE)
		{
			ht->add_kmer(keyF, target, org, position, fstrand);
			cpos--;
		}
	} //it

}

void process_kmergz(gzFile in)
{ 	//read gzip fasta file and extract sequence
	//passes read to process_read
	char buf[BUFLEN];
	char* offset = buf;
	string line, sequence;
        int c1, c2, org, count;
	vtype target;
	ptype position;
	char strand;
	bool fstrand;

	while (true) 
	{
		int err, len = sizeof(buf)-(offset-buf);
		if (len == 0) error("Buffer to small for input line lengths");

 		len = gzread(in, offset, len);

		if (len == 0) break;    
		if (len <  0) error(gzerror(in, &err));

		char* cur = buf;
		char* end = offset+len;

		for (char* eol; (cur<end) && (eol = find(cur, end, '\n')) < end; cur = eol + 1)
 		{
			line = string(cur, eol);
			if (line.back() == '\r')
				line.pop_back();
			if (line.length() > 0)
			{
				replace(line.begin(), line.end(), ',', ' ');
				istringstream ss(line);
				if (ss >> sequence >> target >> org >> position >> strand >> count)
				{
					fstrand = (strand == 'F');
					process_kmer(sequence, target, org, position, fstrand);
					tct++;								
				}
			}
		}
		// any trailing data in [eol, end) now is a partial line
		offset = copy(cur, end, buf);
	}
	// trailing data without eol
	line = string(buf, offset);
	if (gzclose(in) != Z_OK) error("failed gzclose");
	//cout << endl;
}

void process_qual(string acc, string &seq, string &qual)
{ 	//uses PHRED 33 scores
	//to trim read and pass it to process_read
	const int cutoff_qual = 17;
	const char cutoff_char = 32 + cutoff_qual;
	int start, stop, target, i;
	int window_size = 4;	
	int window_val, window_cut = cutoff_qual * window_size;
	string trim_seq;

	stop = seq.length()-1;
	start = 0;
	//trim all low qual bases from end
	while (qual.at(start) < cutoff_char && start < stop)
		start++;
	while (qual.at(stop) < cutoff_char && stop > start)
		stop--;
	//trim by window
	if (start < stop - window_size)
	{
		window_val = 0;
		for (i=0; i<window_size; i++)
			window_val += qual.at(start+i) - 32;
		while (window_val < window_cut && start < stop - window_size)
		{
			window_val += qual.at(start+window_size) - qual.at(start);
			start++;
		}
	}
	if (start < stop - window_size)
	{
		window_val = 0;
		for (i=0; i<window_size; i++)
			window_val += qual.at(stop-i) - 32;
		while (window_val < window_cut && start < stop - window_size)
		{
			window_val += qual.at(stop-window_size) - qual.at(stop);
			stop--;
		}
	}
	//trim_seq = seq.substr(start,stop-start+1);
	if (stop - start >= KSIZE)
	{
		target = process_read(seq,acc,start,stop);
	}

}

void process_fqgz(gzFile in)
{ 	//read gzip fastq file and extract sequence and quality lines
	//passes them to process_qual
	char buf[BUFLEN];
	char* offset = buf;
	string line, seq, acc;
	int mod4=0;

	while (true) 
	{
		int err, len = sizeof(buf)-(offset-buf);
		if (len == 0) error("Buffer to small for input line lengths");

 		len = gzread(in, offset, len);

		if (len == 0) break;    
		if (len <  0) error(gzerror(in, &err));

		char* cur = buf;
		char* end = offset+len;

		for (char* eol; (cur<end) && (eol = find(cur, end, '\n')) < end; cur = eol + 1)
 		{
			line = string(cur, eol);
			if (line.back() == '\r') //windows
				line.pop_back();
			if (line.length() > 0)
			{
				if (mod4 == 1)
				{
					seq = line;
				}
				else if (mod4 == 0)
				{
					acc = line;
				}
				else if (mod4 == 3)
				{
					process_qual(acc, seq, line);
				}
				mod4 = (mod4 + 1) % 4;	
				
			}
			//cout << line << endl;
		}

		// any trailing data in [eol, end) now is a partial line
		offset = copy(cur, end, buf);
	}

	// trailing data without eol
	line = string(buf, offset);

	if (gzclose(in) != Z_OK) error("failed gzclose");
}

void process_fagz(gzFile in)
{ 	//read gzip fasta file and extract sequence
	//passes read to process_read
	char buf[BUFLEN];
	char* offset = buf;
	string line;
        string sequence = "", acc = "";
        int c1, c2, target;

	while (true) 
	{
		int err, len = sizeof(buf)-(offset-buf);
		if (len == 0) error("Buffer to small for input line lengths");

 		len = gzread(in, offset, len);

		if (len == 0) break;    
		if (len <  0) error(gzerror(in, &err));

		char* cur = buf;
		char* end = offset+len;

		for (char* eol; (cur<end) && (eol = find(cur, end, '\n')) < end; cur = eol + 1)
 		{
			line = string(cur, eol);
			if (line.back() == '\r')
				line.pop_back();
			if (line.length() > 0)
			{
				if (line.at(0) == '>')
				{
					if (sequence.length() > KSIZE)
					{ // process previous sequence
						target = process_read(sequence, acc, 0, sequence.length()-1);
					}
					sequence = "";
					acc = line.substr(1,line.length()-1);
				}
				else
				{
					sequence += line;
				}				
			}
			//cout << line << endl;
		}

		// any trailing data in [eol, end) now is a partial line
		offset = copy(cur, end, buf);
	}
	if (sequence.length() > KSIZE)
	{ //last one
		target = process_read(sequence, acc, 0, sequence.length()-1);
	}
	// trailing data without eol
	line = string(buf, offset);
	if (gzclose(in) != Z_OK) error("failed gzclose");

}

void process_fa(string rname)
{	//read unzipped fasta file and extract sequence
	//passes read to process_read
	string sequence = "", line, lseq, acc= "";
	ifstream fin;
	int target;

	fin.open(rname);
	if (fin)
	{
		while (getline(fin, line))
		{
			stringstream linestream(line);
			linestream >> lseq;
			if (lseq[0] == '>')
			{
				if (sequence.length() > KSIZE)
				{
					target = process_read(sequence, acc, 0, sequence.length()-1);
				}
				sequence = "";
				acc = lseq.substr(1,lseq.length()-1);
			}
			else
				sequence += lseq;
		}
		fin.close();
		if (sequence.length() > KSIZE)
		{ //last one
			target = process_read(sequence, acc, 0, sequence.length()-1);
		}
	}
	else
	{
		cout << "nark " << rname << endl;
	}
}

int main(int argc, char* argv[])
{
	int i,j;
	vtype target;
	ptype position;
	char strand;
	bool fstrand;
	string header, fname, line, r1name, r2name, dname, oname1, oname2;
	int refi, targi, strni, count, org;
	ifstream fin;
	string genus, species, acc, sequence, lseq, remainder;
	vector<string> fnames;
	vector<string>::iterator it2;
	string trname;
	string s2 = e1;

	taxonomy = new Tree1();

 	if (argc > 1) 
	{
		dname = argv[1];
 		//if (argc > 2) 
		//	TARGET = atoi(argv[2]);
	}
	else
	{
		dname = "/mnt/dmb/Mark_backup/J/";
	}
	ht = new Hashtable();
	tableM = new int[MAXALIGN*MAXALIGN];
	tableI = new int[MAXALIGN*MAXALIGN];
	tableD = new int[MAXALIGN*MAXALIGN];
	tableinit();
	//load strain list
	for (i = 0; i < MAXTAR; i++)
		ntargorgs[i] = 0;
	fin.open(iname);
	if (fin)
	{
		while (getline(fin, line))
		{
			stringstream linestream(line);
			linestream >> targi >> acc;
			accession[norgs] = acc;
			targno[norgs] = targi;
			if (targi>0)
				ntargorgs[targi]++;
			//cout << norgs << " nt " << targi << " ts " << strni << endl;
			norgs++;
		}
		fin.close();
		//cout << norgs << " strains" << endl;
	}
	else
	{
		cout << "narin " << iname << endl;
	}
	//load taxonomy tree
	fin.open(tname);
	if (fin)
	{
		while (getline(fin, line))
		{
			stringstream linestream(line);
			linestream >> i >> j;
			taxonomy->add_edge(i,j);
		}
	}
	fin.close();
	cout << "tree loaded" << endl;

	//load kmer database
	tct = 0;
	process_kmergz(gzopen(pname.c_str(), "rb"));
	cout << tct << " kmers loaded" << endl;

	//find all matching files in directory
	DIR *dir;
	struct dirent *ent;
	string name1;
	cout << dname << endl;
	if ((dir = opendir (dname.c_str())) != NULL) 
	{
		while ((ent = readdir (dir)) != NULL) 
		{
			name1 = string(ent->d_name);
			size_t pos = name1.find(s2);
			if ( pos != string::npos)
			{
				fnames.push_back(name1.substr(0,pos));
			}
		}
		closedir (dir);
	} else 
	{
		// could not open directory 
		cout << "hosed" << endl;
		perror ("");
		return EXIT_FAILURE;
	}
	for (it2 =fnames.begin(); it2 != fnames.end(); ++it2)
	{
		memset(gcount,0,MAXTAR*sizeof(int));
		memset(ucount,0,MAXTAR*sizeof(int));
                kmer_seen.clear();
		oname2 = dname + *it2 + "_result.txt";
		trname = dname + *it2 + "_reads.txt";
		cout << *it2 << endl;
		tct = 0;
		//if (TARGET > 0)
			outread.open (trname.c_str(), ofstream::out | ofstream::trunc);
#if FASTQ == 1
		r1name = dname + *it2 + e1;
		r2name = dname + *it2 + e2;
		process_fqgz(gzopen(r1name.c_str(), "rb"));
		cout << tct << " reads loaded" << endl;
		process_fqgz(gzopen(r2name.c_str(), "rb"));
#else
		r1name = dname + *it2 + e1;
		process_fa(r1name.c_str());
#endif
		cout << tct << " reads loaded" << endl;

		//if (TARGET > 0)
			outread.close();
		ofstream out2(oname2);
		for (i=0; i<MAXTAR; i++)
			out2 << i << "," << gcount[i] << "," << ucount[i] << endl;
		out2.close();

	}

	delete ht;
	delete[] tableM;
	delete[] tableI;
	delete[] tableD;
	delete taxonomy;

return 0;
}
