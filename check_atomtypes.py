import sys, os

in_file1 = sys.argv[1]
in_file2 = sys.argv[2]

f1 = open(in_file1, 'r')
lines1 = f1.readlines()
f1.close()

f2 = open(in_file2, 'r')
lines2 = f2.readlines()
f2.close()

atypes1 = {}
atypes2 = {}

for line1 in lines1:
    if line1[0] != ";":
        key, anum, mass, charge, ptype, sigma, epsi = line1.split()
        atypes1[key] = {}
        atypes1[key]["anum"] = anum
        atypes1[key]["mass"] = mass
        atypes1[key]["charge"] = charge
        atypes1[key]["ptype"] = ptype
        atypes1[key]["sigma"] = sigma
        atypes1[key]["epsi"] = epsi


for line2 in lines2:
    if line2[0] != ";":
        key, anum, mass, charge, ptype, sigma, epsi = line2.split()
        atypes2[key] = {}
        atypes2[key]["anum"] = anum
        atypes2[key]["mass"] = mass
        atypes2[key]["charge"] = charge
        atypes2[key]["ptype"] = ptype
        atypes2[key]["sigma"] = sigma
        atypes2[key]["epsi"] = epsi

for key1 in atypes1:
    if key1 not in atypes2:
        print(key1 + " is missing from file2")
    print(key1)

for key2 in atypes2:
    if key2 not in atypes1:
        print(key2 + " is missing from file1")
    print(key2)


for key in atypes1:
    if atypes1[key]["anum"] == atypes2[key]["anum"] and atypes1[key]["mass"] == atypes2[key]["mass"] and atypes1[key]["charge"] == atypes2[key]["charge"] and atypes1[key]["ptype"] == atypes2[key]["ptype"] and atypes1[key]["sigma"] == atypes2[key]["sigma"] and atypes1[key]["epsi"] == atypes2[key]["epsi"]:
        print("all good")
    else:
        print(key, atypes1[key], atypes2[key])
