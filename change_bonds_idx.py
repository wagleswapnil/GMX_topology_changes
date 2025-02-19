import sys, os

in_file = sys.argv[1]
idx_file = sys.argv[2]

f = open(in_file, 'r')
lines = f.readlines()
f.close()

fidx = open(idx_file, 'r')
idxs = fidx.readlines()
fidx.close()

idx_to_newidx = {}

for idx in idxs:
    idx_to_newidx[idx.split()[0]] = idx.split()[1]

for line in lines:
    print(idx_to_newidx[line.split()[0]].rjust(7) + "  " + idx_to_newidx[line.split()[1]].rjust(7) + "  " + " ".join(line.split()[2:]))
