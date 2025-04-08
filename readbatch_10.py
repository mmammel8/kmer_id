import re
import sys
from os import listdir
from os.path import isfile, join

dir1 = "/home/mmammel/fastq/"
outname = "test_b10.csv"
ext1 = "_result.txt"

mincount = 2.0
minuniq = 3.0
maxrat = 80000.0

exclude_b = set([])
exclude_i = set([4178,1744,2539,5624,1575,5647,323,2728,268,5317,297,3867,314,1344,2947,2935,4213,4976,2767,2763,118,3390,1757])
exclude_s = set(list(range(1928,2339)))
exclude_b = exclude_b | exclude_i

exclude_s = set(list(range(1928,2339)))
#exclude_s.remove(1931)
exclude_b = exclude_b | exclude_i | exclude_s

reffile = "./bact10/refkey10.txt"
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
#count_arr = np.array(count_list)
num_targs = len(name_list)

#read in all output files in directory
resultfiles = [ f for f in listdir(dir1) if isfile(join(dir1,f)) and f.endswith(ext1)]
num_cols =  len(resultfiles)
noid_list = []
read_ct = []
file_list = []
s_list = []
#m = np.zeros((num_targs,num_cols))
m = [ [0 for _ in range(num_cols)] for _ in range(num_targs)]
col = 0
for f in resultfiles:
    fname = dir1 + f
    pos = f.find(ext1)
    if pos > -1:
    	f = f[:pos]
    file_list.append(f)
    data_file = open(fname, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    read_ct.append(0.0)
    index = 0
    spindex = 0
    for line in lines:
        if len(line) > 1:
            row = line.split(',')
            target = int(row[0])
            count = float(row[1])
            if len(row) > 2:
                uniq = float(row[2])
            else:
                uniq = count
            count2 = count
            if count2 < mincount or uniq < minuniq or (count2 / uniq > maxrat):
                count2 = 0.0
            read_ct[col] += count
            if target > 0:
                if in_use[target] == 1:
                    m[index][col] = count2
                    index += 1
            else:
                noid_list.append(int(count))
    col += 1

#b = m / count_arr[:,None] #normalize each row by kmer count
#sums = np.sum(b, axis=0)
rowmax = [0 for _ in range(num_targs)]
b = [ [0 for _ in range(num_cols)] for _ in range(num_targs)]
sums = [0 for _ in range(num_cols)]
for col in range(num_cols):
    sums[col] = 0
    for row in range(num_targs):
        b[row][col] = m[row][col] / count_list[row]
        sums[col] += b[row][col]
    if sums[col] <  0.00000009:
        sums[col] = 0.0000001
    for row in range(num_targs):
        b[row][col] = b[row][col] * 100.0 / sums[col]
        rowmax[row] = max(rowmax[row], b[row][col])
#b = b / sums[None,:]
#b = b * 100.0
print (num_cols, " x ", num_targs)

order_col = sorted(range(num_cols), key=lambda k: file_list[k])
#order_row = sorted(range(num_targs), key=lambda k: name_list[k])

#rowmax = b.max(axis=1)
#print (rowmax)
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
    if rowmax[i] > 0.000:
        output = name_list[i]
        for j in range(num_cols):
            k = order_col[j]
            output += ',' + str(m[i][k]) + ',' + str(b[i][k])
        output += "\n"
        out_file.write(output)
out_file.close()
