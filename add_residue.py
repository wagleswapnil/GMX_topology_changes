import sys, os, re
import argparse
from protein_top_parser import parse_protein_top
from write_topology import write_protein_topology

def fetch_top_backbone_atoms(protein, idxA):
    nminus1N = None
    nminus1C = None
    nminus1O = None
    nminus1CA = None
    nminus1H = None
    nminus1HA = None
    count = 0
    #This loop picks indices of backbone for n-1 residue, where n is the mutating residue
    for i, atom in enumerate(protein["atoms"]):
        if atom[0] == ";":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if nr == idxA:
                count = i
                break
            if aname ==  "N":
                nminus1N = nr
            elif aname ==  "C":
                nminus1C = nr
            elif aname == "CA":
                nminus1CA = nr
            elif aname == "O":
                nminus1O = nr
            elif aname == "H":
                nminus1H = nr
            elif aname == "HA":
                nminus1HA = nr
            else:
                continue
    print(nminus1N, nminus1C, nminus1H)
    #count = max(int(nminus1N), int(nminus1C), int(nminus1CA), int(nminus1O), int(nminus1H), int(nminus1HA))
    print(count)
    #This loop picks the indices for the nth residue, i.e., the one that is mutating
    #TODO add the case where backbone is also mutating. This loop and the n+1 loop will change
    nN = None
    nC = None
    nCA = None
    nO = None
    nH = None
    nHA = None
    for i, atom in enumerate(protein["atoms"][count:]):
        if atom[0] == ";":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if nN != None and nC != None and nCA != None and nO != None and nH != None and nHA != None:
                count += i
                break
            if aname ==  "N":
                nN = nr
            elif aname ==  "C":
                nC = nr
            elif aname == "CA":
                nCA = nr
            elif aname == "O":
                nO = nr
            elif aname == "H":
                nH = nr
            elif aname == "HA":
                nHA = nr
            else:
                continue
    print(nN, nC, nH)
    #count = max(int(nN), int(nC), int(nCA), int(nO), int(nH), int(nHA))
    print(count)

    #This loop picks indices of backbone for n+1 residue, where n is the mutating residue
    nplus1N = None
    nplus1C = None
    nplus1CA = None
    nplus1O = None
    nplus1H = None
    nplus1HA = None
    for i, atom in enumerate(protein["atoms"][count:]):
        if atom[0] == ";":
            continue
        else:
            if nplus1N != None and nplus1C != None and nplus1CA != None and nplus1O != None and nplus1H != None and nplus1HA != None:
                break
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if aname ==  "N":
                nplus1N = nr
            elif aname ==  "C":
                nplus1C = nr
            elif aname == "CA":
                nplus1CA = nr
            elif aname == "O":
                nplus1O = nr
            elif aname == "H":
                nplus1H = nr
            elif aname == "HA":
                nplus1HA = nr
            else:
                continue
    print(nplus1N, nplus1C, nplus1H)
    return

def add_residue_connections(proteinA, idxA, proteinB):
    beginA = None
    endA = None
    bb_atoms = ["N", "C", "O", "CA", "H", "HA", "CB"]
    aname_to_idx = {}
    if "-" in idxA:
        beginA = int(idxA.split("-")[0])
        endA = int(idxA.split("-")[1])
    else:
        beginA = int(idxA.strip())
        endA = int(idxA.strip())

    for i, atom in enumerate(proteinA["atoms"]):
        if atom[0] == ";" or atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) < beginA:
                if aname in bb_atoms:
                    key = "nminus1" + aname
                    aname_to_idx[key] = nr
            elif beginA <= int(nr) <= endA:
                key = "n_" + residue + "_" + aname
                aname_to_idx[key] = nr
            else:
                break
    for i, atom in enumerate(reversed(proteinA["atoms"])):
        if atom[0] == ";" or atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) > endA:
                if aname in bb_atoms:
                    key = "nplus1" + aname
                    aname_to_idx[key] = nr
            else:
                break
    
    proteinA["bonds"].append("; bonds for the dummy residue")
    for bond in proteinB["bonds"]:
        ai, aj, rest = bond.split(None, 2)
        proteinA["bonds"].append(aname_to_idx[ai] + "    " + aname_to_idx[aj] + "    " + rest)

    return aname_to_idx, proteinA

def get_residue_connections(proteinB, idxB):
    beginB = None
    endB = None
    bb_atoms = ["N", "C", "O", "CA", "H", "HA", "CB"]
    real_bb_atoms = ["N", "C", "O", "CA", "H", "HA"]
    real_bb_params = {}
    idx_to_aname = {}
    if "-" in idxB:
        beginB = int(idxB.split("-")[0])
        endB = int(idxB.split("-")[1])
    else:
        beginB = int(idxB.strip())
        endB = int(idxB.strip())
    
    for i, atom in enumerate(proteinB["atoms"]):
        if atom[0] == ";":
            continue
        elif atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) < beginB:
                if aname in bb_atoms:
                    key = "nminus1" + aname
                    idx_to_aname[nr] =  key
            elif int(nr) >= beginB and int(nr) <= endB:
                key = "n_" + residue + "_" + aname
                idx_to_aname[nr] =  key
            else:
                break
    for i, atom in enumerate(reversed(proteinB["atoms"])):
        if atom[0] == ";":
            continue
        elif atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) > endB:
                if aname in bb_atoms:
                    key = "nplus1" + aname
                    idx_to_aname[nr] =  key
            else:
                break
    r = range(beginB, endB+1)

    bonds = []
    for bond in proteinB["bonds"]:
        if bond[0] == ";" or bond[0] == "#":
            continue
        ai, aj, rest = bond.split(None, 2)
        if int(ai) in r or int(aj) in r:
            bonds.append(idx_to_aname[ai] + "    " + idx_to_aname[aj] + "    " + rest)
    proteinB["bonds"] = bonds

    pairs = []
    for pair in proteinB["pairs"]:
        if pair[0] == ";" or pair[0] == "#":
            continue
        ai, aj, rest = pair.split(None, 2)
        if int(ai) in r or int(aj) in r:
            pairs.append(idx_to_aname[ai] + "    " + idx_to_aname[aj] + "    " + rest)
    proteinB["pairs"] = pairs

    angles = []
    for angle in proteinB["angles"]:
        if angle[0] == ";" or angle[0] == "#":
            continue
        ai, aj, ak, rest = angle.split(None, 3)
        if int(ai) in r or int(aj) in r or int(ak) in r:
           angles.append(idx_to_aname[ai] + "    " + idx_to_aname[aj] + "    " + idx_to_aname[ak] + "    " + rest) 
    proteinB["angles"] = angles

    dihedrals = []
    for dihedral in proteinB["dihedrals"]:
        if dihedral[0] == ";" or dihedral[0] == "#":
            continue
        ai, aj, ak, al, rest = dihedral.split(None, 4)
        if int(ai) in r or int(aj) in r or int(ak) in r or int(al) in r:
            dihedrals.append(idx_to_aname[ai] + "    " + idx_to_aname[aj] + "    " + idx_to_aname[ak] + "    " + idx_to_aname[al] + "    " + rest)
    proteinB["dihedrals"] = dihedrals
    return idx_to_aname, proteinB

def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--topA", type=str, help="path for topology file A, i.e., the WT topology", required=True)
    parser.add_argument("--topB", type=str, help="path for topology file A, i.e., the mutant topology", required=True)
    parser.add_argument("--idxA", type=str, help="index of atom right before the mutating residue", required=True)
    parser.add_argument("--idxB", type=str, help="range of indices to pick from mutant to add to the WT", required=True)
    parser.add_argument("--backbone", type=str, help="should the backbone be added as well, yes/no?", required=False)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypesB, proteinB = parse_protein_top(args.topB, "system1")
    atomtypesA, proteinA = parse_protein_top(args.topA, "system1")
    idx_to_aname, proteinB = get_residue_connections(proteinB, args.idxB)
    aname_to_idx, proteinA = add_residue_connections(proteinA, args.idxA, proteinB)
    for key, value in aname_to_idx.items():
        print(key + "\t" + value)
    print(proteinA["bonds"])
    return

if __name__== "__main__":
    main()
    
