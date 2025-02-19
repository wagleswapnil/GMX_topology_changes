import sys, os

in_file = sys.argv[1]

f = open(in_file, 'r')
lines = f.readlines()
f.close()

ridx = 224
idx = 3572

for i in range(1, idx + 1):
    print(str(i) + "  " +str(i))

for line in lines:
    if "residue" in line:
        ridx += 1
#        print("; residue  " + str(ridx) + " " + " ".join(line.split()[3:]))
    else:
        idx += 1
#        print(" " + str(idx) + "         " + line.split()[1] + "    " +  str(ridx) + "    " + line.split()[3] + " " + line.split()[4].rjust(6) + "    " + "  ".join(line.split()[5:]))
        print(line.split()[0] + "  " + str(idx))
