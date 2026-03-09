#Reads kmer outputs in a directory and generates report
#MKM 5/10/2016

import re
import sys
from os import listdir
from os.path import isfile, join
import numpy as np

mypath = "W:/Mark_backup/ROAR/Saffron/chloro/"
outname ="saffron_chloro.csv"
mincount = 2.0
minuniq = 2.0
maxrat = 80.0

exclude_b = []

reffile = "refKeyc3.txt"
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
        target,name,count,hit,tested, gsize, nstrains  = line.split('\t')
        row = name.split('_')
        target = int(target)
        hit = float(hit)
        count = float(count)
        if nstrains != '0':
            gensize = float(gsize) / float(nstrains)
        else:
            gensize = 1.0
        tested = float(tested)
        use = "1"
        if target in exclude_b or count < 10.0 or hit < 10.0 or len(row) < 6:
            use = "0"
        in_use.append(int(use))
        if use == "1":
            name_list.append(name)
            factor_list.append(tested / hit / gensize)
factor_arr = np.array(factor_list)
num_targs = len(name_list)

#read in all output files in directory
onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) and f[-11:] == "_result.txt" ]
num_cols =  len(onlyfiles)
noid_list = []
read_ct = []
file_list = []
m = np.zeros((num_targs,num_cols))
u = np.zeros((num_targs,num_cols))
col = 0
for f in onlyfiles:
    fname = mypath + '/' + f
    sname = f[:-11]
    file_list.append(sname)
    data_file = open(fname, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    read_ct.append(0.0)
    index = 0
    for line in lines:
        if len(line) > 1:
            row = line.split(',')
            target = int(row[0])
            count = float(row[1])
            uniq = float(row[2])
            count2 = count
            if count2 < mincount or uniq < minuniq or (count2 / uniq > maxrat):
                count2 = 0.0
            read_ct[col] += float(count)
            if target > 0:
                if in_use[target] == 1:
                    m[index, col] = float(count2)
                    u[index, col] = uniq
                    index += 1
            else:
                noid_list.append(int(count))
    col += 1

b = m * factor_arr[:,None] #normalize each row by kmer coverage
sums = np.sum(b, axis=0)
for col in range(num_cols):
    if sums[col] <  0.00000009:
        sums[col] = 0.0000001
print (sums)
b = b / sums[None,:]
b = b * 100.0
print (num_cols, num_targs)

order_col = sorted(range(num_cols), key=lambda k: file_list[k])
#order_row = sorted(range(num_targs), key=lambda k: name_list[k])

rowmax = b.max(axis=1)
uniqmax = u.max(axis=1)
print (rowmax)
out_file = open(outname, 'w')
output = "name,"
for i in range(num_cols):
    output += file_list[order_col[i]] + ",,"
output += "\n"
out_file.write(output)
output = "total,"
for i in range(num_cols):
    output += str(read_ct[order_col[i]]) + ",,"
output += "\n"
out_file.write(output)
output = "no_id,"
for i in range(num_cols):
    output += str(noid_list[order_col[i]]) + ",,"
output += "\n"
out_file.write(output)
for i in range(num_targs):
    #l = order_row[i]
    if rowmax[i] > 0.000:
        output = name_list[i]
        for j in range(num_cols):
            k = order_col[j]
            #output += ',' + str(m[l,k]) + ',' + str(b[l,k])
            output += ',' + str(m[i,k]) + ',' + str(b[i,k])
        output += "\n"
        out_file.write(output)
out_file.close()
