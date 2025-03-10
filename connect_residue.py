import sys, os, re
import argparse
from protein_top_parser import parse_protein_top
from write_topology import write_protein_topology
from get_b_connections import get_resB_connections, get_topB_indices
from insert_b_connections import get_topA_indices, make_resB_connections

def get_topB_connections(proteinB, resB):
    topB_bb_indices, resB_indices, resnameB = get_topB_indices(proteinB, resB)
    proteinB = get_resB_connections(proteinB, topB_bb_indices, resB_indices)
    return proteinB, resnameB

def insert_stateB_connections(proteinA, resA, stateB_connections, resnameB):
    topA_indices, resA_indices = get_topA_indices(proteinA, resA, resnameB)
    for key, value in topA_indices.items():
        print(key + "\t" + value)
    proteinA = make_resB_connections(proteinA, stateB_connections, topA_indices, resA_indices)
    return proteinA


def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--topA", type=str, help="path for topology file A, i.e., the WT topology", required=True)
    parser.add_argument("--topB", type=str, help="path for topology file A, i.e., the mutant topology", required=True)
    parser.add_argument("--resA", type=str, help="index of mutating residue", required=True)
    parser.add_argument("--resB", type=str, help="index of residue to pick from mutant to add to the WT", required=False)
    parser.add_argument("--out-path", type=str, help="filename for output", required=False)
    parser.add_argument("--backbone", type=str, help="should the backbone be added as well, yes/no?", required=False)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypesB, proteinB = parse_protein_top(args.topB, "system1")
    atomtypesA, proteinA = parse_protein_top(args.topA, "system1")
    if args.resB == None:
        args.resB = args.resA
    stateB_connections, resnameB = get_topB_connections(proteinB, args.resB)
    print(resnameB)
    #print("".join(stateB_connections["dihedrals"]))
    proteinA = insert_stateB_connections(proteinA, args.resA, stateB_connections, resnameB)
    write_protein_topology(atomtypesA, proteinA, args.out_path)
    return

if __name__== "__main__":
    main()

