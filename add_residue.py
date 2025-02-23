import sys, os, re
import argparse
from protein_top_parser import parse_protein_top

def fetch_topA_atoms(protein, idx):
    nminus1N = None
    for i, atom in enumerate(protein["atoms"]):
        if atom[0] == ";":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if 



    return


def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--topA", type=str, help="path for topology file A, i.e., the WT topology", required=True)
    parser.add_argument("--topB", type=str, help="path for topology file A, i.e., the mutant topology", required=False)
    parser.add_argument("--idxA", type=str, help="index of atom after which new residue is to be added", required=False)
    parser.add_argument("--idxB", type=str, help="range of indices to pick from mutant to add to the WT", required=False)
    parser.add_argument("--backbone", type=str, help="should the backbone be added as well, yes/no?", required=False)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypes, protein = parse_protein_top(args.top, "system1")
    if args.A:
        atomtypes, protein = generate_Bstate(atomtypes, protein, args.A)
    if args.B:
        atomtypes, protein = generate_Astate(atomtypes, protein, args.B)
    return

if __name__== "__main__":
    main()
    
