import sys, os, re
import argparse
from protein_top_parser import parse_protein_top

def write_protein_topology(atomtypes, protein):
    moleculetype_section = ["moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals", "cmap", "virtual_sites2", "constraints", "exclusions", "other", "settles"]
    print("[ atomtypes ]")
    print("".join(atomtypes))
    for key in moleculetype_section:
        if key in protein:
            print(f"[ {key} ]")
            print("".join(protein[key]))
    return
def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--top", type=str, help="path for the topology file", required=True)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypes, protein = parse_protein_top(args.top, "system1")
    write_protein_topology(atomtypes, protein)
    return

if __name__== "__main__":
    main()

