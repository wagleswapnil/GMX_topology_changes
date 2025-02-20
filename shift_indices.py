import sys, os, re
import argparse
from protein_top_parser import parse_protein_top

def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--topA", type=str, help="path for topology A, i.e., the topology of the wild type protein", required=True)
    parser.add_argument("--topB", type=str, help="path for topology B, i.e., the topology of the mutant protein", required=True)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    print(args)
    parse_protein_top(args.topA, "system1")
    return

if __name__== "__main__":
    main()
    
