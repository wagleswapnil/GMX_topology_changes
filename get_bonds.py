import sys, os

in_file = sys.argv[1]
idx_file = sys.argv[2]
uidx_file = sys.argv[3]

f = open(in_file, 'r')
lines = f.readlines()
f.close()

fidx = open(idx_file, 'r')
idxs = fidx.readlines()
fidx.close()

uidx = open(uidx_file, 'r')
uidxs = uidx.readlines()
uidx.close()

idx_list = []
for idx in idxs:
    if idx.strip() in idx_list:
        pass
    else:
        idx_list.append(idx.strip())


uidx_dict = {}
for uidx in uidxs:
    uidx_dict[uidx.split()[0]] = uidx.split()[4]

for line in lines:
    if line.split()[0] in idx_list or line.split()[1] in idx_list or line.split()[2] in idx_list or line.split()[3] in idx_list:
        if line.split()[2] == "1":
            print(uidx_dict[line.split()[0]] + "    " + uidx_dict[line.split()[1]] + " " + " ".join(line.split()[2:]))
        elif line.split()[3] == "1":
            #print(line.split())
            print(uidx_dict[line.split()[0]] + "    " + uidx_dict[line.split()[1]] + "    " + uidx_dict[line.split()[2]] + " " + " ".join(line.split()[3:]))
        elif line.split()[4] == "1":
            print(uidx_dict[line.split()[0]] + "    " + uidx_dict[line.split()[1]] + "    " + uidx_dict[line.split()[2]] + "    " + uidx_dict[line.split()[3]] + " " + " ".join(line.split()[4:]))
