import sys, os, re
import argparse
from protein_top_parser import parse_protein_top
from write_topology import write_protein_topology

def get_atomtypes_dict(atomtypes):
    atomtypes_dict = {}
    for i, atomtype in enumerate(atomtypes):
        if atomtype[0] == ";":
            continue
        else:
            atomtypes_dict[atomtype.split()[0]] = atomtype.split(None,1)[1].strip()
    return atomtypes_dict

def generate_Bstate(atomtypes, protein, A, backbone):
    atomtypes_dict = get_atomtypes_dict(atomtypes)
    bb_atoms = ["N", "C", "O", "CA", "H", "HA"]
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
                if backbone == "no":
                    if aname in bb_atoms:
                        continue
                protein["atoms"][i] = nr.rjust(6) + "  " + atype.rjust(6) + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "  " + aname.rjust(6) + "  " + cgnr.rjust(6) + "  " + q.rjust(12) + "  " + m.rjust(10) + "    " + "DUM" + atype + "    " + "0.000000" + "    " + m.rjust(10) + " " + rest
                if "DUM" + atype in atomtypes_dict:
                    continue
                else:
                    atomtypes_dict["DUM" + atype] = "0.000000     0.000000   A     0.000000     0.000000"
                    atomtypes.append("DUM" + atype + "    " + "0.000000     0.000000   A     0.000000     0.000000\n")
    return atomtypes, protein

def generate_Astate(atomtypes, protein, B, backbone):
    atomtypes_dict = get_atomtypes_dict(atomtypes)
    bb_atoms = ["N", "C", "O", "CA", "H", "HA"]
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
                if backbone == "no":
                    if aname in bb_atoms:
                        continue
                protein["atoms"][i] = nr.rjust(6) + "   " + "DUM" + atype + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "    " + "D" + aname + "  " + cgnr.rjust(6) + "  " + "0.000000" + "    " + m.rjust(10) + "    " + atype.rjust(6) + "    " + q.rjust(12) + "  " + m.rjust(10) + " " + rest
                if "DUM" + atype in atomtypes_dict:
                    continue
                else:
                    atomtypes_dict["DUM" + atype] = "0.000000     0.000000   A     0.000000     0.000000"
                    atomtypes.append("DUM" + atype + "    " + "0.000000     0.000000   A     0.000000     0.000000\n")
    return atomtypes, protein

def check_atomtypes():
    return

def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--topA", type=str, help="path for topology file A", required=True)
    parser.add_argument("--topB", type=str, help="path for topology file B", required=True)
    parser.add_argument("--A", type=str, help="range of indices, for which A state is interacting and B state is dummy", required=True)
    parser.add_argument("--B", type=str, help="range of indices, for which B state is interacting and A state is dummy", required=True)
    parser.add_argument("--backbone", type=str, help="should A and B states  be generated for the backbone, yes/no?", required=True)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypesA, proteinA = parse_protein_top(args.topA, "system1")
    atomtypesB, proteinB = parse_protein_top(args.topB, "system1")
    atomtypes = check_atomtypes(atomtypesA, atomtypesB)
    protein = merge_atoms(proteinA, args.A, proteinB, args.B)
    if args.A:
        atomtypes, protein = generate_Bstate(atomtypes, protein, args.A, args.backbone)
    if args.B:
        atomtypes, protein = generate_Astate(atomtypes, protein, args.B, args.backbone)
    write_protein_topology(atomtypes, protein)
    return

if __name__== "__main__":
    main()
    
