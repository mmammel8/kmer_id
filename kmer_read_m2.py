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
    foldername = wdir + "/mitochondria"
    in2 = foldername + "_data.txt"
    in3 = foldername + "_key.txt"
    in4 = foldername + "_count.txt"
    cname = wdir + "/" + "result.txt"
    wdir += "/"

    #read in reference names and counts
    count_list = []
    name_list = []
    name_dict = dict()
    in_use = []

    data_file = open(in3, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    for line in lines:
        if len(line) > 1:
            target,name = line.split('\t')
            name_dict[target] = name
     
    data_file = open(in4, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    for line in lines:
        if len(line) > 1:
            target,count = line.split(',')
            use =  (int(count) > 35)
            in_use.append(use)
            if use:
                name_list.append(name_dict[target])
                count_list.append(float(count)+10.0) #not too low

    count_arr = np.array(count_list)
    num_targs = len(name_list)

    process = Popen([wdir + '/kmerread', '-wdir', wdir, '-f1', file1, '-f2', file2])
    (stdout, stderr) = process.communicate()

    #read in output files
    noid_list = []
    num_cols = 1
    m = np.zeros((num_targs,num_cols))
    col = 0
    data_file = open(cname, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    index = 0
    for line in lines:
        if len(line) > 1:
            t_s, count, uniq = line.split(',')
            target = int(t_s)
            if target > 0:
                if in_use[target]: 
                    m[index, col] = float(count)
                    index += 1
            else:
                noid_list.append(int(count))

    b = m / count_arr[:,None] #normalize each row by kmer count
    sums = np.sum(b, axis=0)
    b = b / sums[None,:]
    b = b * 100.0
    rowmax = b.max(axis=1)

    out_file = open(oname, 'w')
    output = "taxid,reads,abundance\n"
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




