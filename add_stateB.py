import sys, os, re
import argparse
from protein_top_parser import parse_protein_top

def get_atomtypes_dict(atomtypes):
    atomtypes_dict = {}
    for i, atomtype in enumerate(atomtypes):
        if atomtype[0] == ";":
            continue
        else:
            atomtypes_dict[atomtype.split()[0]] = atomtype.split(None,1)[1].strip()
    return atomtypes_dict

def generate_Bstate(atomtypes, protein, A):
    atomtypes_dict = get_atomtypes_dict(atomtypes)
    beginA = None
    endA = None
    if "-" in A:
        beginA = int(A.split("-")[0])
        endA = int(A.split("-")[1])
    else:
        beginA = int(A.strip())
        endA = int(A.strip())
    for i, atom in enumerate(protein["atoms"]):
        if atom[0] == ";":
            continue
        elif atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) >= beginA and int(nr) <= endA:
                protein["atoms"][i] = nr.rjust(6) + "  " + atype.rjust(6) + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "  " + aname.rjust(6) + "  " + cgnr.rjust(6) + "  " + q.rjust(12) + "  " + m.rjust(10) + "    " + "DUM" + atype + "    " + "0.000000" + "    " + m.rjust(10) + " " + rest
                if "DUM" + atype in atomtypes_dict:
                    continue
                else:
                    atomtypes_dict["DUM" + atype] = "0.000000     0.000000   A     0.000000     0.000000"
                    atomtypes.append("DUM" + atype + "    " + "0.000000     0.000000   A     0.000000     0.000000\n")
    return atomtypes, protein

def generate_Astate(atomtypes, protein, B):
    atomtypes_dict = get_atomtypes_dict(atomtypes)
    beginB = None
    endB = None
    if "-" in B:
        beginB = int(B.split("-")[0])
        endB = int(B.split("-")[1])
    else:
        beginB = int(B.strip())
        endB = int(B.strip())
    for i, atom in enumerate(protein["atoms"]):
        if atom[0] == ";":
            continue
        elif atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) >= beginB and int(nr) <= endB:
                protein["atoms"][i] = nr.rjust(6) + "   " + "DUM" + atype + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "    " + "D" + aname + "  " + cgnr.rjust(6) + "  " + "0.000000" + "    " + m.rjust(10) + "    " + atype.rjust(6) + "    " + q.rjust(12) + "  " + m.rjust(10) + " " + rest
                if "DUM" + atype in atomtypes_dict:
                    continue
                else:
                    atomtypes_dict["DUM" + atype] = "0.000000     0.000000   A     0.000000     0.000000"
                    atomtypes.append("DUM" + atype + "    " + "0.000000     0.000000   A     0.000000     0.000000\n")
    return atomtypes, protein

def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--top", type=str, help="path for topology file", required=True)
    parser.add_argument("--A", type=str, help="range of indices, for which A state is interacting and B state is dummy", required=False)
    parser.add_argument("--B", type=str, help="range of indices, for which B state is interacting and A state is dummy", required=False)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypes, protein = parse_protein_top(args.top, "system1")
    if args.A:
        atomtypes, protein = generate_Bstate(atomtypes, protein, args.A)
    if args.B:
        atomtypes, protein = generate_Astate(atomtypes, protein, args.B)
    print("".join(atomtypes))
    print("".join(protein["atoms"][10:20]))
    return

if __name__== "__main__":
    main()
    
