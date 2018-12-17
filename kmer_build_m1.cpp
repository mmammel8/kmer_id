// newkmer_3.cpp find k-mers for all species
// save in one file
// MKM 1 APR 2016
//-std=c++0x -D__NO_INLINE__
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
#include <time.h>
#include <sys/time.h>
#include <zlib.h>
#include <cstddef>
#include <cstring>
#include <algorithm>
#include <set>
using namespace std;

//data files
string iname = "data.txt";   //text file of groups for probe design
string tname = "tree.txt";   //edges in taxonomy tree
string iname2 = "filter.txt";   //text file of groups for exclusion
string oname = "probes.txt";

typedef unsigned long long ui64;
typedef unsigned long long ktype;
typedef unsigned int vtype;
typedef unsigned long long itype;

const int KSIZE = 30; //kmer size
//const int MAXORGS = 10855;
//const int MAXTAR = 17136; //2 to MAXTAR-1 inclusive are present, 1 is root, 0 unused
const itype MAXHASH = (1LL << 34); //34 //pow2
const int MAXREP = 2048;
const int REPSHIFT = 11; //(log of maxrep)
const int TARGET = 0;
const int MAXPROBES = 80000;
const int MINTARG = 0; //for pruning

const char bases[4] = { 'A', 'C', 'G', 'T' };
const ui64 one = 1LL;
const ui64 two = 2LL;
const ui64 three = 3LL;
const ui64 mask = (one << (KSIZE * 2)) - 1;
const ui64 hiC = one << (KSIZE - 1) * 2;
const ui64 hiG = two << (KSIZE - 1) * 2;
const ui64 hiT = three << (KSIZE - 1) * 2;
const ui64 hi1 = one << 62;
const ui64 hi2 = two << 62;
const ui64 hi3 = three << 62;
const ui64 himask = hi1 - one;
vector<int> ntargorgs, pcount;
int mcount, tct;

string int2str(ktype key)
{
	string sequence = "";
	for (int i = 0; i < KSIZE; ++i)
	{
		sequence = bases[key & 3] + sequence;
		key >>= 2;
	}
	return sequence;
}

class Tree1
{
public:

	//vector<int> children[MAXTAR];
	vector<int> parent;

	Tree1(int nt)
	{
		for (int i = 0; i < nt; i++)
			parent.push_back(1); //root
	}

	~Tree1()
	{

	}

	void add_edge(int x, int y)
	{
		if (x < parent.size() && y < parent.size())
			parent[y] = x;
		//children[x].push_back(y);
	}

	int ca(int x, int y)
	{
		//return most recent common ancestor of nodes x and y
		set<int> ancestors;
		int z;

		ancestors.insert(1);
		z = x;
		while (z > 1)
		{
			ancestors.insert(z);
			z = parent[z];
		}
		z = y;
		while (ancestors.find(z) == ancestors.end())
		{
			z = parent[z];
		}
		return z;
	}

	int get_parent(int x)
	{
		if (x>1)
			return parent[x];
		else
			return 1;
	}

};

Tree1 *taxonomy;

class Hashtable
{
public:

	itype size;
	vtype *pHash;

	Hashtable()
	{
		size = 0;
		pHash = new vtype[MAXHASH];
		HashClear();
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
		memset(pHash, 0, MAXHASH * sizeof(vtype));
	}

	void HashAdd(ktype key, int targi)
	{ 
		itype index;
		int target, count;

		index = (integerHash(key)) & (MAXHASH-1); //for power of 2 size
		target = (pHash[index] >> REPSHIFT);
		if (pHash[index] == 0)
		{ //empty cell
			//pHash[index] = targi;
			pHash[index] = ((targi << REPSHIFT) | 1);
			if (++size > MAXHASH - 32)
			{
				cout << "out of memory in table " << endl;
				exit(1);
			}
		}
		else if (target > 1)
		{
			target = taxonomy->ca(target, targi);
			count = (pHash[index] & (MAXREP-1));
			if (count == (MAXREP - 1))
			{
				//cout << "too many reps " << target << endl;
				pHash[index] = 1; //0 target = bad
			}
			else
				pHash[index] = ((target << REPSHIFT) | (count+1));
		}
	}

	void HashRemove(ktype key)
	{ 
		itype index;

		index = (integerHash(key)) & (MAXHASH-1); //for power of 2 size
		if (pHash[index] > 1)
		{ 
			pHash[index] = 1; //0 target = bad
		}
	}

	vtype getHash(ktype key, int &count)
	{
		itype index = (integerHash(key)) & (MAXHASH - 1);
		vtype target = (pHash[index] >> REPSHIFT);
		count = (pHash[index] & (MAXREP-1));
		pHash[index] = 1; //DO NOT REUSE
		return target;
	}

}; //class Hashtable
Hashtable *ht;

string process_fa(string seq)
{  
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

const unsigned BUFLEN = 0x4000;

void error(const char* const msg)
{
    cerr << msg << "\n";
    exit(255);
}

string process_gz(gzFile in)
{ //read gzip fastq file and extract sequence
	char buf[BUFLEN];
	char* offset = buf;
	string line;
        string sequenc = "";

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
					sequenc += "N"; //contig separator
				}
				else
				{
					sequenc += process_fa(line);
				}				
			}
		}
		// any trailing data in [eol, end) now is a partial line
		offset = copy(cur, end, buf);
	}

	// trailing data without eol
	line = string(buf, offset);
	if (gzclose(in) != Z_OK) error("failed gzclose");

	return sequenc;
}

void process_seq(string sequence, int targi)
{ //add family kmer to tree
	char base;
	string::iterator it1;
	int cpos, mpos;
	ktype keyF, keyR;

	cpos = 0;
	mpos = 0;
	keyF = 0;
	keyR = 0;
	for (it1 = sequence.begin(); it1 != sequence.end(); ++it1)
	{
		base = *it1;
        mpos++;
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
		default:
			cpos = 0; //ambiguous character
			mpos = 0;
			keyF = 0;
			keyR = 0;
			break;
		}
		if (cpos == KSIZE)
		{
			ht->HashAdd(keyF < keyR ? keyF : keyR, targi);
			cpos--;
		}
	} //it

}

void process_seq3(string sequence, int targi)
{ //remove kmers matching outside of family
	char base;
	string::iterator it1;
	int cpos, mpos;
	ktype keyF, keyR;

	cpos = 0;
	mpos = 0;
	keyF = 0;
	keyR = 0;
	for (it1 = sequence.begin(); it1 != sequence.end(); ++it1)
	{
		base = *it1;
        mpos++;
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
		default:
			cpos = 0; //ambiguous character
			mpos = 0;
			keyF = 0;
			keyR = 0;
			break;
		}
		if (cpos == KSIZE)
		{
			ht->HashRemove(keyF < keyR ? keyF : keyR);
			cpos--;
		}
	} //it

}


bool check_entropy(string sequence)
{
	int i, row, maxrow;
	double as, cs, gs, ts, total, pa, pc, pg, pt, entropy;
	char base, prev = 'N';
	string::iterator it1;

	as = 1.0;
	cs = 1.0;
	gs = 1.0;
	ts = 1.0;
	row = 0;
	maxrow = 0;
	for (it1 = sequence.begin(); it1 != sequence.end(); ++it1)
	{
		base = *it1;
		if (prev == base)
		{
			if (++row > maxrow)
				maxrow = row;
		}
		else
		{
			row = 1;
			prev = base;
		}
		switch (base)
		{
		case 'A':
			as+=1.0;
			break;
		case 'C':
			cs+= 1.0;
			break;
		case 'G':
			gs+= 1.0;
			break;
		case 'T':
			ts+= 1.0;
			break;
		}
	}
	if (maxrow > 11) //for KSIZE = 25/30
		return false;
	total = as + cs + gs + ts;
	pa = as/total;
	pc = cs/total;
	pg = gs/total;
	pt = ts/total;
	entropy = -pa*log10(pa) - pc*log10(pc) - pg*log10(pg) - pt*log10(pt);
	//cout << entropy << endl;
	if (entropy < 0.43) //for KSIZE = 25/30
		return false;

	return true;

}

int process_seq2(string sequence, int org, int targi, ofstream &out)
{
	//save kmers to file
	char base;
	string::iterator it1;
	int cpos, mpos, gpos, position, minpos=-1;
	ktype keyF, keyR, key;
	int target, count, minct, nprobes=0;
	string kmer;
	char fstrand;

	cpos = 0;
	mpos = 0;
	keyF = 0;
	keyR = 0;
	gpos = 0;
	for (it1 = sequence.begin(); it1 != sequence.end(); ++it1)
	{
		base = *it1;
		mpos++;
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
		default:
			cpos = 0; //ambiguous character
			mpos = 0;
			keyF = 0;
			keyR = 0;
			break;
		}
		if (cpos == KSIZE)
		{
			key = keyF < keyR ? keyF : keyR;
			target = ht->getHash(key, count);//DOES NOT REUSE
			if (keyF < keyR)
				fstrand = 'F';
			else
				fstrand = 'R';
			position = gpos - KSIZE + 1; //start position
			if (ntargorgs[target] < 3)
			    minct = ntargorgs[target];
			else if (ntargorgs[target] < 6)
				minct = ntargorgs[target] - 1;
			else
				minct = ntargorgs[target] / 5 + 1;
			if (target < MINTARG) minct = 1;
			if ((target > 1) && (count >= minct) && (mpos > minpos) && (pcount[target] < MAXPROBES))
			{
				kmer = int2str(key);
				if (check_entropy(kmer))
				{
					out << kmer << "," << target << "," << org << "," << gpos << "," << fstrand << "," << count << endl;
					if (targi >= MINTARG)
						minpos = mpos + KSIZE; //non-overlapping
					else
						minpos = mpos + 2;
					pcount[target]++;
					tct++;
				}
			}
			cpos--;
		}
		gpos++;
	} //it
	return 0;

}

int main(int argc, char* argv[])
{
	int i, j, n;
	vector<string> accession, accession2;
	vector<int> targno;
	string header, fname, line, sequence;
	int num_orgs=0, num_orgs2=0, targi, num_targ=0;
	string argst, fdir = "", wdir = "", name = "bob";

 	if (argc > 1) 
	{
		for (i = 1; i < argc; i++)
		{
			argst = argv[i];
    			if (argst == "-name")
			{
				name = argv[i+1];
				wdir = "./" + name + "/";
			}
    			if (argst == "-fadir")
			{
				fdir = argv[i+1];
			}
 		}
	}
	iname = wdir + name + "_" + iname; 
	tname = wdir + name + "_" + tname;  
	iname2 = wdir + name + "_" + iname2; 
	oname = wdir + name + "_" + oname;

	ifstream fin;
	string acc;
	ofstream out(oname.c_str());
	ht = new Hashtable();

	//load access of outgroup sequences
	fin.open(iname2);
	if (fin)
	{
		while (getline(fin, line))
		{
			stringstream linestream(line);
			if (line.back() == '\r')
				line.pop_back();
			linestream >> acc;
			accession2.push_back(acc);
			num_orgs2++;
		}
	}
	fin.close();
	cout << num_orgs2 << " outs loaded" << endl;
	cout.flush();
	//load ingroups
	num_targ = 0;
	fin.open(iname);
	if (fin)
	{
		while (getline(fin, line))
		{
			stringstream linestream(line);
			if (line.back() == '\r')
				line.pop_back();
			linestream >> targi >> acc;
			accession.push_back(acc);
			targno.push_back(targi);
			if (targi > num_targ)
				num_targ = targi;
			num_orgs++;
		}
	}
	fin.close();
	num_targ++;
	cout << num_orgs << " sequences loaded" << endl;
	ntargorgs.assign(num_targ,0);
	pcount.assign(num_targ,0);
	taxonomy = new Tree1(num_targ);
	for (i = 0; i < num_orgs; i++)
	{ //count number of strains for each target
		targi = targno[i];
		while (targi > 1)
		{
			ntargorgs[targi]++;
			targi = taxonomy->get_parent(targi);
		}
	}

	//load tree
	fin.open(tname);
	if (fin)
	{
		while (getline(fin, line))
		{
			stringstream linestream(line);
			if (line.back() == '\r')
				line.pop_back();
			linestream >> i >> j;
			taxonomy->add_edge(i,j);
		}
	}
	fin.close();
	cout << "tree loaded" << endl;

	for (i = 0; i < num_orgs; i++)
	{
		if (targno[i] > 1)
		{
			fname = fdir + accession[i] + ".fasta.gz";
			cout << "1 " << i << " " << num_orgs << " " << accession[i] << "\r";
			sequence = process_gz(gzopen(fname.c_str(), "rb"));
			process_seq(sequence, targno[i]);
		}
	}
	cout << endl;
	//outgroups
	for (i = 0; i < num_orgs2; i++)
	{
		fname = fdir + accession2[i] + ".fasta.gz";
		cout << "2 " << i << " " << num_orgs2 << " " << accession2[i] << "\r";
		sequence = process_gz(gzopen(fname.c_str(), "rb"));
		process_seq3(sequence, 1);

	}
	cout << endl;
	//save good ones
	tct = 0;
	for (i = 0; i < num_orgs; i++)
		if (targno[i] > 1)
		{
			fname = fdir + accession[i] + ".fasta.gz";
			cout << "3 " << i << " " << num_orgs << " " << accession[i] << "\r";
			sequence = process_gz(gzopen(fname.c_str(), "rb"));
			sequence = process_fa(sequence);
			process_seq2(sequence, i, targno[i], out);
		}
	out.close();
	cout << endl;
	oname = wdir + name + "_count.txt";
	ofstream out2(oname.c_str());
	for (i=0; i<num_targ; i++)
	{
		out2 << i << "," << pcount[i] << endl;
	}
	out2.close();
	cout << "probe count " << tct << endl;
	cout << "size " << ht->size << endl;

	delete ht;
	delete taxonomy;

return 0;
}
