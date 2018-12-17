#Create taxonomy tree
#MKM 10/10/2016
#root is 1, 0 is unused

import sys
from subprocess import Popen, PIPE

infile = "mitochondria_list.txt"
num_cols = 6 #taxonomy columns
foldername = "mitochondria"
fastadir = "/mnt/dmb/Mark_backup/mitochondria/"

dir = "./" + foldername + "/"
outname = dir + foldername + "_tree.txt"
out2 = dir + foldername + "_data.txt"
out3 = dir + foldername + "_key.txt"
infile = dir + infile
name_dict = dict()
tree = dict()
data_file = open(infile, 'r')
data = ''.join(data_file.readlines())
data_file.close()
lines = data.split('\n')

next_tax = 2
curr_tax = 0
for line in lines:
    row = line.split('\t')
    if len(row) > num_cols and row[0] != "class":
        prev_tax = 0
        for pos in range(num_cols):
            if len(row[pos]) > 0:
                name = "_".join(row[:pos+1])
                if name in name_dict:
                    curr_tax = name_dict[name]
                else:
                    curr_tax = next_tax
                    name_dict[name] = next_tax
                    next_tax += 1
                if prev_tax > 0:
                    if prev_tax in tree:
                        tree[prev_tax].add(curr_tax)
                    else:
                        tree[prev_tax] = set([curr_tax])
                prev_tax = curr_tax
            
out_file = open(outname, 'w')
for prev_tax in tree:
    for curr_tax in tree[prev_tax]:
        output = str(prev_tax) + "\t" + str(curr_tax) + "\n"
        if prev_tax == curr_tax:
            print prev_tax, "!!!!"
        out_file.write(output)
out_file.close()

out_file = open(out2, 'w')
for line in lines:
    row = line.split('\t')
    if len(row) > num_cols and row[0] != "class":
        curr_tax = 0
        for pos in range(num_cols):
            if len(row[pos]) > 0:
                name = "_".join(row[:pos+1])
                curr_tax = name_dict[name]
        output = str(curr_tax) + "\t" + row[num_cols] + "\n"
        out_file.write(output)
out_file.close()

out_file = open(out3, 'w')
for word in name_dict:
    curr_tax = name_dict[word]
    output = str(curr_tax) + "\t" + word + "\n"
    out_file.write(output)
out_file.close()

process = Popen(['./kmerbuild', '-name', foldername, '-fadir', fastadir])
(stdout, stderr) = process.communicate()

