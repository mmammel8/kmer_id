// newkmer_11.cpp load k-mers for all species
// test reads by SmithWaterman
// MKM 1 APR 2018
//-std=c++0x ... at end: -lz
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

#define FAST 0
#define FASTQ 1  //fasta = 0
unordered_set<int> myset = { 1928,1929,1930,1931,1932,2048,2051,2107,2178,2309,2330 };

typedef int ptype;
typedef unsigned long long ui64;
typedef unsigned long long ktype;
typedef unsigned int vtype;
typedef unsigned short int otype;
typedef unsigned long long itype;

int TARGET = 0; //save these, second input parameter
const int KSIZE = 30; //kmer size
const int MAXORGS = 14791;
const int MAXTAR = 5982; //2 to MAXTAR-1 inclusive are present, 1 is root, 0 unused
const int MAXREP = 2048;
const int REPSHIFT = 11; //(log of maxrep)
const itype MAXHASH = (1 << 30); //pow2
//const itype MAXHASH = 1048576*64;
const int GAPO =11;
const int GAPX =1;
const int MATCH = 5;
const int MISMATCH = -4;
const bool INIGAPPEN = false;
//const int threshold_scr = 900;
string accession[MAXORGS];
int targno[MAXORGS], ntargorgs[MAXTAR], mcount;
int gcount[MAXTAR], savepos, tct, norgs=0;
ofstream outread;

//input data files
string iname = "bData10.txt";   //text file of groups for probe design
string pname = "./probes7/probes10.txt";
string tname = "btree_10.txt";   //edges in taxonomy tree
string dir = "/mnt/dmb/Mark_backup/genbank/";

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
		while (z != root)
		{
			ancestors.insert(z);
			z = get_parent(z);
		}
		z = y;
		if (ancestors.find(y) != ancestors.end())
			return x; //x is descendent of y
		while (ancestors.find(z) == ancestors.end())
		{
			z = get_parent(z);
			if (z==x) //y is descendent of x
				return y;
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

string load_data3(string fnameIn3)
{
	//loads sequence from FASTA file, concats contigs and separates them with 'N'
	ifstream fin;
	string sequence, seq2;
	string::iterator it1;
	char base;

	seq2 = "";
	fin.open(fnameIn3);
	if (!fin.is_open())
	{
		cout << "could not find " << fnameIn3 << endl;
		return seq2;
	}
	while (getline(fin, sequence))
	{
		sequence.erase(remove_if(sequence.begin(), sequence.end(), ::isspace), sequence.end());
		if (sequence.length() > 1)
		{
			if (sequence.at(0) == '>')
			{
				seq2 += "N"; //contig separator
			}
			else
			{
				for (it1 = sequence.begin(); it1 != sequence.end(); ++it1)
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
			}
		}
	}
	fin.close();

	return seq2;
}

int align(string dna1, string dna2)
{
	int *tableM, *tableI, *tableD;
	int b,bb,c,cc, len1, len2;
	char a1,a2;
	int i, j, max, rn, startx, starty;
	const int NINF = -2000000000;

	len1 = dna1.length() + 1;
	len2 = dna2.length() + 1;
	tableM = new int[len1*len2];
	tableI = new int[len1*len2];
	tableD = new int[len1*len2];
	tableM[0]=0;
	tableI[0]=NINF;
	tableD[0]=NINF;
	for (i=1; i<len1; i++) //initialize upper row
	{
		tableI[0+i]=NINF;
		if (INIGAPPEN) 
		{
			tableD[0+i] = -GAPO - (i-1)*GAPX;
			tableM[0+i] = -GAPO - (i-1)*GAPX;
		}
		else 
		{
			tableD[0+i] = 0;
			tableM[0+i] = 0;
		}
	}
	for (j=1; j<len2; j++) //init left column
	{
		tableD[j*len1+0]=NINF;
		if (INIGAPPEN) 
		{
			tableI[j*len1+0] = -GAPO - (j-1)*GAPX;
			tableM[j*len1+0] = -GAPO - (j-1)*GAPX;
		}
		else 
		{
			tableI[j*len1+0] = 0;
			tableM[j*len1+0] = 0;
		}
	}
	j=1;
	do
	{
		i=1;
		do
		{
			a1=dna1[i-1];
			a2=dna2[j-1];
			max=tableM[(j-1)*len1 + i-1];
			if (tableI[(j-1)*len1 + i-1] > max)
				max=tableI[(j-1)*len1 + i-1];
			if (tableD[(j-1)*len1 + i-1] > max)
				max=tableD[(j-1)*len1 + i-1];
			//copy/replace
                        if (a1 != a2)
                            max += MISMATCH;
			else max += MATCH;
			tableM[j*len1 + i] = max;

			b=tableM[(j-1)*len1 + i] - GAPO; //new gap
			bb = tableI[(j-1)*len1 + i] -GAPX; //extend
			b = (bb<b) ? b : bb;
			tableI[j*len1 + i] = b;

			c=tableM[j*len1 + i-1] - GAPO; //new gap
			cc = tableD[j*len1 + i-1] -GAPX; //extend
			c = (cc<c) ? c : cc;
			tableD[j*len1 + i] = c;		
		} while(++i < len1);
	} while(++j < len2);

	startx=len1-1;
	starty=len2-1;
	max=tableM[starty*len1 + startx]; //last entry;
	if (tableI[starty*len1 + startx] > max)
		max=tableI[starty*len1 + startx];
	if (tableD[starty*len1 + startx] > max)
		max=tableD[starty*len1 + startx];
	delete[] tableM;
	delete[] tableI;
	delete[] tableD;

	//printf("cost = %d\n",max);

	return max;

}

int process_read(string &sequence, int start, int stop)
{
	//returns id of read, 0 = no match
	char base;
	int it1, it2;
	string kmer, sequence2, dna1, dna1r, dna2, fname;
	ktype keyF, keyR, key;
	int i, cpos, maxscr;
	bool fstrand1, fstrand2;
	int readlength = stop - start + 1;
	int st2, scr, stlen2;
	ktype target, prev_targ=0;
	ptype position;
	otype org;
	int maxct, totct;

	cpos = 0;
	keyF = 0;
	keyR = 0;
	maxscr = 5 * sequence.length() / 2;
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

#if FAST==0
			if (target > 0 && gcount[target] < 100 && target != prev_targ)
			{ //do alignment on initial reads of a species
				fname = dir + accession[org] + ".fasta";
				sequence2 = load_data3(fname);
				stlen2 = sequence2.length();
				fstrand1 = (keyF < keyR);
				if (fstrand1 == fstrand2)
				{
					st2 = position + KSIZE - it1 - 1;
					if (st2 < 0) st2 = 0;
					if (st2 + readlength >= stlen2) st2 = stlen2 - readlength;
					dna1 = sequence.substr(start, readlength);
					dna2 = sequence2.substr(st2, readlength);
				}
				else
				{
					st2 = position + it1 - readlength + 1;
					if (st2 < 0) st2 = 0;
					if (st2 + readlength >= stlen2) st2 = stlen2 - readlength;
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
					dna2 = sequence2.substr(st2, readlength);
					//reverse complement dna1
				}
				scr = align(dna1, dna2);
				if (scr < maxscr)
				{
					target = 1; //fail
				}
			}
#endif
			if (prev_targ > 0)
			{
				target = taxonomy->msca(target, prev_targ);
			}
			prev_targ = target;
			cpos--;
		}
	} //it

	return target;
}

void process_kmer(string sequence, vtype target, otype org, ptype position, bool fstrand)
{
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

void process_qual(string acc, string seq, string qual)
{ //uses PHRED 33 scores
	const int cutoff_qual = 17;
	const char cutoff_char = 32 + cutoff_qual;
	int start, stop, target;
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
		window_val = qual.at(start) + qual.at(start+1) + qual.at(start+2) + qual.at(start+3);
		while (window_val < window_cut && start < stop - window_size)
		{
			window_val += qual.at(start+4) - qual.at(start);
			start++;
		}
	}
	if (start < stop - window_size)
	{
		window_val = qual.at(stop) + qual.at(stop-1) + qual.at(stop-2) + qual.at(stop-3);
		while (window_val < window_cut && start < stop - window_size)
		{
			window_val += qual.at(stop-4) - qual.at(stop);
			stop--;
		}
	}
	//trim_seq = seq.substr(start,stop-start+1);
	if (stop - start >= KSIZE)
	{
		target = process_read(seq,start,stop);
		//unordered_set<int>::const_iterator got = myset.find(target);
		//if ( got != myset.end() )
		if (target == TARGET && TARGET>0)
		{
			//string s = to_string(savepos);
	 		outread << acc << endl << seq.substr(start, stop-start+1) << endl;
		}
		gcount[target]++;
		tct++; 
	}

}

const unsigned BUFLEN = 0x4000;

void error(const char* const msg)
{
    cerr << msg << "\n";
    exit(255);
}

void process_gz(gzFile in)
{ //read gzip fastq file and extract sequence and quality lines
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

void process_fa(string rname)
{
	string sequence = "", line, lseq;
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
					process_read(sequence, 0, sequence.length()-1);
					tct++;
				}
				sequence = "";
			}
			else
				sequence += lseq;
		}
		fin.close();
		if (sequence.length() > KSIZE)
		{ //last one
			target = process_read(sequence, 0, sequence.length()-1);
			gcount[target]++;
			tct++;
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
#if FASTQ == 1
	string s2 = "_L001_R1_001.fastq.gz";
#else
	string s2 = "_R1.fasta";
#endif
	taxonomy = new Tree1();

 	if (argc > 1) 
	{
		dname = argv[1];
 		if (argc > 2) 
			TARGET = atoi(argv[2]);
	}
	else
	{
		dname = "/mnt/dmb/Susan/LeafyGreens/Cilantro/PreEnrichment/";
	}
	ht = new Hashtable();

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
		cout << norgs << " sequences loaded" << endl;
	}
	else
	{
		cout << "narin " << iname << endl;
	}
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

	tct = 0;
	ifstream infile(pname);
	while (infile)
	{
		string s1;
		if (!getline(infile, s1)) break;
		replace(s1.begin(), s1.end(), ',', ' ');
		istringstream ss(s1);
		if (ss >> sequence >> target >> org >> position >> strand >> count)
		{
			fstrand = (strand == 'F');
			process_kmer(sequence, target, org, position, fstrand);
			tct++;
		}
	}
	cout << tct << " kmers loaded" << endl;

/*
sequence = "TCTGGTGGCCAGCCTTAATGATGTTCGGGCCAAACGATGATAACTCGCCAAACAGCGCCAGAAGTCTCACCTGGAAAATCAAACGTTTCACCAACGACGAACTCCGCCAGCGTTTCGTGGATAACACCGTTCCA";
process_read(sequence,0,133);
sequence = "GCCACAACTTGTTGATGGCCTGCTGCATCTTTTGCCCTGATACGTCAGTACCATTGCCCAGTCGCTCCAGCCAGCCACGACTAAAACGCAGGTGATAGCGCGCTTCTTTAATTGCTTTGGCAGAAATCGCCGCCAGTTGCGGATCACGGCTTTCCATCAGACGGGTAAAGAGC";
process_read(sequence,0,172);
*/

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
		oname2 = dname + *it2 + "_result.txt";
		trname = dname + *it2 + "_reads.txt";
		cout << *it2 << endl;
		tct = 0;
		if (TARGET > 0)
			outread.open (trname.c_str(), ofstream::out | ofstream::trunc);
#if FASTQ == 1
		r1name = dname + *it2 + "_L001_R1_001.fastq.gz";
		r2name = dname + *it2 + "_L001_R2_001.fastq.gz";
		process_gz(gzopen(r1name.c_str(), "rb"));
		cout << tct << " reads loaded" << endl;
		process_gz(gzopen(r2name.c_str(), "rb"));
#else
		r1name = dname + *it2 + "_R1.fasta";
		process_fa(r1name.c_str());
#endif
		cout << tct << " reads loaded" << endl;

		if (TARGET > 0)
			outread.close();
		ofstream out2(oname2);
		for (i=0; i<MAXTAR; i++)
			out2 << i << "," << gcount[i] << endl;
		out2.close();


	}

	delete ht;

return 0;
}
