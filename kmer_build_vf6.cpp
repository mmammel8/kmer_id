// newkmer_3.cpp find k-mers for all species
// MKM 1 APR 2016
//-std=c++0x -D__NO_INLINE__

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
#include <sys/stat.h>
#include <unistd.h>
using namespace std;

//data files
string iname = "data.txt";   //text file of groups for probe design
string tname = "tree.txt";   //edges in taxonomy tree
string iname2 = "filter.txt";   //text file of groups for exclusion
string oname = "probes.txt";

typedef unsigned long long ui64;
typedef unsigned long long ktype; //key
typedef unsigned int vtype; //value

const int KSIZE = 30; //kmer size must be no more than 1/2 bits in ktype
const ktype MAXHASH = (1LL << 35); //34 //pow2
const int MAXREP = 2048; //leaves 2,000,000 possible taxids. 2048 is too small?
const int REPSHIFT = 11; //(log of maxrep)
const int TARGET = 0;
const int MAXPROBES = 100000;
const int MINTARG = 0; //for pruning
//string outdir = "/mnt/dmb/Mark_backup/genbank_fung/ncbi-genomes/"; //make variable?
//string outdir = "/home/MMammel/newkmer/mitochondria5/genbank_fung/"; //make variable?
string outdir = "/mnt/dmb/Mark_backup/genbank/";


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

	ktype size;
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
	ktype integerHash(ktype k)
	{
		k ^= k >> 33;
		k *= 0xff51afd7ed558ccd;
		k ^= k >> 33;
		k *= 0xc4ceb9fe1a85ec53;
		k ^= k >> 33;
		return (ktype) k;
	}

	void HashClear()
	{
		memset(pHash, 0, MAXHASH * sizeof(vtype));
	}

	void HashAdd(ktype key, int targi)
	{ 
		ktype index;
		int target, count;

		index = (integerHash(key)) & (MAXHASH-1); //for power of 2 size
		target = (pHash[index] >> REPSHIFT);
		if (pHash[index] == 0)
		{ //empty cell
			//pHash[index] = targi;
			pHash[index] = ((targi << REPSHIFT) | 1);
			++size;
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
		ktype index;

		index = (integerHash(key)) & (MAXHASH-1); //for power of 2 size
		if (pHash[index] > 1)
		{ 
			pHash[index] = 1; //0 target = bad
		}
	}

	vtype getHash(ktype key, int &count)
	{
		ktype index = (integerHash(key)) & (MAXHASH - 1);
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

string load_data2(string fnameIn3)
{
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
				seq2 += process_fa(sequence);
			}
		}
	}
	fin.close();

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
	double total[10], pa[10], pc[10], pg[10], pt[10], entropy[10];
	char base, prev = 'N';
	string::iterator it1;
	double log10_4 = log10(4.0);
	bool result;

	ktype keyF, keyR = 0;

	for (i=0; i<10; i++)
	{
		pa[i] = 1.0;
		pc[i] = 1.0;
		pg[i] = 1.0;
		pt[i] = 1.0;
		total[i] = 0.0;
	}
	row = 0;
	maxrow = 0;
	i = 0;
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
			keyF = (keyF << 2) & mask;
 			pa[i%2] += 1.0;
			pa[i%3+2] += 1.0;
			pa[i%5+5] += 1.0;
			break;
		case 'C':
			keyF = ((keyF << 2) & mask) | 1;
 			pc[i%2] += 1.0;
			pc[i%3+2] += 1.0;
			pc[i%5+5] += 1.0;
			break;
		case 'G':
			keyF = ((keyF << 2) & mask) | 2;
 			pg[i%2] += 1.0;
			pg[i%3+2] += 1.0;
			pg[i%5+5] += 1.0;
			break;
		case 'T':
			keyF = ((keyF << 2) & mask) | 3;
 			pt[i%2] += 1.0;
			pt[i%3+2] += 1.0;
			pt[i%5+5] += 1.0;
			break;
		}
		i++;
	}
	if (maxrow > 11) //for KSIZE = 25/30
		return false;
	for (i=0; i<10; i++)
		total[i] = pa[i] + pc[i] + pg[i] + pt[i];
	for (i=0; i<10; i++)
	{
		pa[i] = pa[i]/total[i];
		pc[i] = pc[i]/total[i];
		pg[i] = pg[i]/total[i];
		pt[i] = pt[i]/total[i];
		entropy[i] = -pa[i]*log10(pa[i]) - pc[i]*log10(pc[i]) - pg[i]*log10(pg[i]) - pt[i]*log10(pt[i]);
	}

	double e2 = (entropy[0] + entropy[1]) / 2.0 /log10_4;
	double e3 = (entropy[2] + entropy[3] + entropy[4]) / 3.0 /log10_4;
	double e5 = (entropy[5] + entropy[6] + entropy[7]  + entropy[8] + entropy[9]) / 5.0 /log10_4;

	//cout << e2 << "\t" << e3 << "\t" << e5 << endl;

	if (e2 < 0.80 || e3 < 0.80 || e5 < 0.80)
		return false;

  	bool bad = ((keyF & 0x3333333333333333) == 0LL) || ((keyF & 0xCCCCCCCCCCCCCCCC) == 0LL);
	if (bad)
		cout << sequence << endl;

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
			if (ntargorgs[target] == 1)
 				minct = 1;
			else if (ntargorgs[target] < 4)
				minct = 2;
			else if (ntargorgs[target] < 10)
				minct = ntargorgs[target] - 2;
			else
				minct = ntargorgs[target] / 5 + 1;
			//if (target < MINTARG) minct = 1;
			if ((target > 1) && (count >= minct) && (gpos > minpos) && (pcount[target] < MAXPROBES))
			{
				kmer = int2str(key);
				if (check_entropy(kmer))
				{
					out << kmer << "," << target << "," << org << "," << gpos << "," << fstrand << "," << count << endl;
					//if (targi >= MINTARG)
					minpos = gpos + KSIZE; //non-overlapping
					//else
					//minpos = mpos + 1;
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

inline bool file_exists (const string& name) 
{
	struct stat buffer;   
	return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc, char* argv[])
{
	int i, j, n;
	vector<string> accession, accession2;
	vector<int> targno;
	string header, fname,fname1,fname2,fname3, line, sequence;
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
			fname1 = fdir + accession[i] + ".fasta.gz";
			fname2 = outdir + accession[i] + ".fasta.gz";
			fname3 = fdir + accession[i] + "_contigs.fasta";
			if (file_exists(fname1))
			{
				sequence = process_gz(gzopen(fname1.c_str(), "rb"));
			}
			else if (file_exists(fname2))
			{
				sequence = process_gz(gzopen(fname2.c_str(), "rb"));
			}
			else if (file_exists(fname3))
			{
				sequence = load_data2(fname3);
			}
			else
			{
				cout << "no file for " << accession[i] << endl;
				exit(1);
			}
			cout << "1 " << i << " " << num_orgs << " " << accession[i] << "\n";
			process_seq(sequence, targno[i]);
		}
	}
	cout << endl;

	//outgroups
	for (i = 0; i < num_orgs2; i++)
	{
		//cout << accession2[i] << endl;
		fname1 = outdir + accession2[i] + ".fasta.gz";
		fname2 = fdir + accession2[i] + ".fasta.gz";
		if (file_exists(fname1))
		{
			sequence = process_gz(gzopen(fname1.c_str(), "rb"));
		}
		else if (file_exists(fname2))
		{
			sequence = process_gz(gzopen(fname2.c_str(), "rb"));
		}
		else
		{
			cout << "no file for " << accession2[i] << endl;
			exit(1);
		}
		cout << "2 " << i << " " << num_orgs2 << " " << accession2[i] << "\n";
		process_seq3(sequence, 1);

	}
	cout << endl;
	//save good ones
	tct = 0;
	for (i = 0; i < num_orgs; i++)
		if (targno[i] > 1)
		{
			fname1 = fdir + accession[i] + ".fasta.gz";
			fname2 = outdir + accession[i] + ".fna.gz";
			fname3 = fdir + accession[i] + "_contigs.fasta";
			if (file_exists(fname1))
			{
				sequence = process_gz(gzopen(fname1.c_str(), "rb"));
			}
			else if (file_exists(fname2))
			{
				sequence = process_gz(gzopen(fname2.c_str(), "rb"));
			}
			else if (file_exists(fname3))
			{
				sequence = load_data2(fname3);
			}
			else
			{
				cout << "no file for " << accession[i] << endl;
				exit(1);
			}
			cout << "3 " << i << " " << num_orgs << " " << accession[i] << "\n";
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
