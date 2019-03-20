#Reads fastq filenames, runs kmerread, and generates report
#MKM 1/1/2019

import re
import sys
import numpy as np
from subprocess import Popen, PIPE

if __name__ == '__main__':
    if len(sys.argv) > 2:
        item = sys.argv[1].strip()
        dirnm = sys.argv[2].strip()
        if item == "-w":
            wdir = dirnm
        else:
            print "w?"
    if len(sys.argv) > 4:
        item = sys.argv[3].strip()
        dirnm = sys.argv[4].strip()
        if item == "-d":
            oname = dirnm + "/mitokmer_result.csv"
        else:
            print "d?"

    if len(sys.argv) > 7:
        item = sys.argv[5].strip()
        filenm1 = sys.argv[6].strip()
        filenm2 = sys.argv[7].strip()
        if item == "-i":
            file1 = filenm1
            file2 = filenm2
        else:
            print "i?"

    target = 0
    wdir += "/"
    reffile = wdir + "mitochondria_refkey.txt"
    cname = wdir + "result.txt"
    exclude_b = []

    #read in reference names and counts
    factor_list = []
    name_list = []
    in_use = []
    data_file = open(reffile, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    lines.pop(0) #header
    for line in lines:
        if len(line) > 1:
            target, name, count, hit, tested, gsize, nstrains  = line.split('\t')
            row = name.split('_')
            target = int(target)
            hit = float(hit)
            if nstrains != '0':
	        gensize = float(gsize) / float(nstrains)
            else:
                gensize = 1.0
            tested = float(tested)
            use = "1"
            if target in exclude_b or count < 10.0 or hit < 10.0 or len(row) < 5:
                use = "0"
            in_use.append(int(use))
            if use == "1":
                name_list.append(name)
                factor_list.append(tested / hit / gensize)
    factor_arr = np.array(factor_list)
    num_targs = len(name_list)
    process = Popen([wdir + '/kmerread', '-wdir', wdir, '-f1', file1, '-f2', file2])

    (stdout, stderr) = process.communicate()

    #read in output files
    noid_list = []
    read_ct = []
    num_cols = 1
    m = np.zeros((num_targs,num_cols))
    col = 0
    data_file = open(cname, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    read_ct.append(0.0)
    index = 0
    for line in lines:
        if len(line) > 1:
            t_s, count, uniq = line.split(',')
            target = int(t_s)
            read_ct[col] += float(count) 
            if target > 0:
                if in_use[target]: 
                    m[index, col] = float(count)
                    index += 1
            else:
                noid_list.append(int(count))

    b = m * factor_arr[:,None] #normalize each row by kmer coverage
    sums = np.sum(b, axis=0)
    b = b / sums[None,:]
    b = b * 100.0
    rowmax = b.max(axis=1)

    out_file = open(oname, 'w')
    output = "taxid,reads,abundance\n"
    out_file.write(output)
    output = "total,"
    for i in range(num_cols):
        output += str(read_ct[i]) + ",,"
    output += "\n"
    out_file.write(output)
    output = "no_id,"
    for i in range(num_cols):
        output += str(noid_list[i]) + ",,"
    output += "\n"
    out_file.write(output)
    for i in range(num_targs):
        #l = order_row[i]
        if rowmax[i] > 0.000: #only output non-zero results
            output = name_list[i]
            for j in range(num_cols):
                output += ',' + "{0:.0f}".format(m[i,j]) + ',' + "{0:.3f}".format(b[i,j])
            output += "\n"
            out_file.write(output)
    out_file.close()




