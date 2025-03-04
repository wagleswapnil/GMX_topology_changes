import sys, os, re
import argparse
from protein_top_parser import parse_protein_top
from write_topology import write_protein_topology


def get_topB_indices(proteinB, resB):
    residues = {str(int(resB) - 1): "nminus1", resB: "n", str(int(resB) + 1): "nplus1"}
    bb_atoms = ["N", "CA", "C", "O", "H", "HA"]
    topB_bb_indices = {}
    for atom in proteinB["atoms"]:
        if atom[0] == ";" or atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if resnr in residues and aname in bb_atoms:
                topB_bb_indices[nr] = residues[resnr] + aname
            elif resnr == resB:
                topB_bb_indices[nr] = residues[resnr] + aname
    return topB_bb_indices


def get_topB_connections(proteinA, proteinB, resA, resB):
    topB_bb_indices = get_topB_bb_indices(proteinB, resB)
    for key, value in topB_bb_indices.items():
        print(key + "\t" + value)
    return


def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--topA", type=str, help="path for topology file A, i.e., the WT topology", required=True)
    parser.add_argument("--topB", type=str, help="path for topology file A, i.e., the mutant topology", required=True)
    parser.add_argument("--resA", type=str, help="index of mutating residue", required=True)
    parser.add_argument("--resB", type=str, help="index of residue to pick from mutant to add to the WT", required=False)
    parser.add_argument("--backbone", type=str, help="should the backbone be added as well, yes/no?", required=False)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypesB, proteinB = parse_protein_top(args.topB, "system1")
    atomtypesA, proteinA = parse_protein_top(args.topA, "system1")
    if args.resB == None:
        args.resB = args.resA
    stateB_connections = get_topB_connections(proteinA, proteinB, args.resA, args.resB)
    quit()
    real_bb_params, idx_to_aname, proteinB = get_residue_connections(proteinB, args.resB)
    print(real_bb_params)
    for key, value in idx_to_aname.items():
        print(key + "\t" + value)
    quit()
    aname_to_idx, proteinA = add_residue_connections(proteinA, args.idxA, proteinB, real_bb_params)
    for key, value in aname_to_idx.items():
        print(key + "\t" + value)
    print(proteinA["bonds"])
    return

if __name__== "__main__":
    main()

