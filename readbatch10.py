#Reads kmer outputs in a directory and generates report
#MKM 5/10/2016

import re
import sys
from os import listdir
from os.path import isfile, join
import numpy as np

#mypath = "W:/Susan/Coconut"
mypath = "W:/Susan"
#mypath = "W:/TammyMiSeq/Carmen/FloraUdoLot3_BestLot1"
outname = "TimeZero_b.csv"

exclude_b = [] #[56,3722,5601]

reffile = "refkey10.txt"
#read in reference names and counts
count_list = []
name_list = []
in_use = []
data_file = open(reffile, 'r')
data = ''.join(data_file.readlines())
data_file.close()
lines = data.split('\n')
lines.pop(0) #header
for line in lines:
    if len(line) > 1:
        target,name,count,use  = line.split('\t')
        if int(target) in exclude_b:
            use = "0"
        #if count == "0":
        #    count = "1" #don't divide by 0
        in_use.append(int(use))
        if use == "1":
            count_list.append(float(count)+10.0) #not too low
            name_list.append(name)
count_arr = np.array(count_list)
num_targs = len(name_list)

#read in all output files in directory
onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) and f[-11:] == "_result.txt" ]
num_cols =  len(onlyfiles)
noid_list = []
file_list = []
m = np.zeros((num_targs,num_cols))
col = 0
for f in onlyfiles:
    fname = mypath + '/' + f
    sname = f[:-11]
    file_list.append(sname)
    data_file = open(fname, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    index = 0
    for line in lines:
        if len(line) > 1:
            t_s, count = line.split(',')
            target = int(t_s)
            if target > 0:
                if in_use[target] == 1: 
                    m[index, col] = float(count)
                    index += 1
            else:
                noid_list.append(int(count))
    col += 1

b = m / count_arr[:,None] #normalize each row by kmer count
sums = np.sum(b, axis=0)
b = b / sums[None,:]
b = b * 100.0
print num_cols, num_targs, 

order_col = sorted(range(num_cols), key=lambda k: file_list[k])
#order_row = sorted(range(num_targs), key=lambda k: name_list[k])

rowmax = b.max(axis=1)
print rowmax
out_file = open(outname, 'w')
output = "name,"
for i in range(num_cols):
    output += file_list[order_col[i]] + ",,"
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

"""
if __name__ == '__main__':
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
    else:
        print 'This test requires an input file.  Please select one from the data directory. '
"""