#!/usr/bin/python
"""
Description: This script writes the topology of a molecule in a file. If an input file is given, it reads the topology of moleculetype
"system1" from the input file and writes it in the output file.

Written on: March 06th, 2025
Written by: Swapnil Wagle, swapnil.wagle92[at]gmail[.]com
"""
import sys, os, re
import argparse
from protein_top_parser import parse_protein_top

#Function to write the topology fle
def write_protein_topology(atomtypes, protein, out_top):
    moleculetype_sections = ["moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals", "cmap", "virtual_sites2", "constraints", "exclusions", "other", "settles"]
    with open(out_top, 'w') as f:
        out_top.write("[ atomtypes ]" + "\n")
        out_top.write("".join(atomtypes))
    for key in moleculetype_sections:
        if key in protein:
            out_top.write(f"[ {key} ]\n")
            out_top.write("".join(protein[key]))
        f.close()
    return

#Input argument parser
def input_data():
    parser = argparse.ArgumentParser(description="This script write a GROMACS topology of a molecule in a file.")
    parser.add_argument("--top", type=str, help="path for the input topology file", required=True)
    parser.add_argument("--out-top", type=str, help="path for the output topology file", required=False, default="out.top")
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypes, protein = parse_protein_top(args.top, "system1")
    write_protein_topology(atomtypes, protein, args.out_top)
    return

if __name__== "__main__":
    main()

