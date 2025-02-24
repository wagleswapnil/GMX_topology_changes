import sys, os, re
import argparse
from protein_top_parser import parse_protein_top

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


def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--topA", type=str, help="path for topology file A, i.e., the WT topology", required=True)
    parser.add_argument("--topB", type=str, help="path for topology file A, i.e., the mutant topology", required=False)
    parser.add_argument("--idxA", type=str, help="index of atom right before the mutating residue", required=False)
    parser.add_argument("--idxB", type=str, help="range of indices to pick from mutant to add to the WT", required=False)
    parser.add_argument("--backbone", type=str, help="should the backbone be added as well, yes/no?", required=False)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypes, protein = parse_protein_top(args.topA, "system1")
    fetch_top_backbone_atoms(protein, args.idxA)
    fetch_connections(args.topB, )
    return

if __name__== "__main__":
    main()
    
