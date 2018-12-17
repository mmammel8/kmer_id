#Reads kmer outputs in a directory and generates report
#MKM 2/12/2018

import re
import sys
import numpy as np
from subprocess import Popen, PIPE

jobs_name = "jobs3m"
target = 0
foldername = "mitochondria"
fastadir = "/mnt/dmb/Mark_backup/mitochondria/"

dird = "./" + foldername + "/"
in2 = dird + foldername + "_data.txt"
in3 = dird + foldername + "_key.txt"
in4 = dird + foldername + "_count.txt"
in1 =     "./" + jobs_name + "/" + jobs_name + ".txt"
outname = "./" + jobs_name + "/" + jobs_name + ".csv"

#read in reference names and counts
count_list = []
name_list = []
name_dict = dict()
in_use = []
job_list = []

data_file = open(in3, 'r')
data = ''.join(data_file.readlines())
data_file.close()
lines = data.split('\n')
for line in lines:
    if len(line) > 1:
        target,name = line.split('\t')
        name_dict[target] = name

data_file = open(in1, 'r')
data = ''.join(data_file.readlines())
data_file.close()
lines = data.split('\n')
skip = 0
for line in lines:
    if len(line) > 1:
        if skip == 0:
            jname, skip = line.split()
            job_list.append(jname)
            skip = int(skip)
        else:
            skip -= 1
 
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

process = Popen(['./kmerread', '-name', foldername, '-fadir', fastadir, '-jname', jobs_name, '-target', target])
(stdout, stderr) = process.communicate()

#read in all output files
num_cols =  len(job_list)
noid_list = []
m = np.zeros((num_targs,num_cols))
col = 0
for f in job_list:
    fname = "./" + jobs_name + "/" + f + "_result.txt"
    data_file = open(fname, 'r')
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
    col += 1

b = m / count_arr[:,None] #normalize each row by kmer count
sums = np.sum(b, axis=0)
b = b / sums[None,:]
b = b * 100.0

order_col = sorted(range(num_cols), key=lambda k: job_list[k])
#order_row = sorted(range(num_targs), key=lambda k: name_list[k])

rowmax = b.max(axis=1)
#print rowmax
out_file = open(outname, 'w')
output = "name,"
for i in range(num_cols):
    output += job_list[order_col[i]] + ",,"
output += "\n"
out_file.write(output)
output = "no_id,"
for i in range(num_cols):
    output += str(noid_list[order_col[i]]) + ",,"
output += "\n"
out_file.write(output)
for i in range(num_targs):
    #l = order_row[i]
    if rowmax[i] > 0.000: #only output non-zero results
        output = name_list[i]
        for j in range(num_cols):
            k = order_col[j]
            output += ',' + str(m[i,k]) + ',' + str(b[i,k])
        output += "\n"
        out_file.write(output)
out_file.close()

"""
if __name__ == '__main__':
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
    else:
        print 'This test requires an input file.  Please select one from the data directory. '
"""
