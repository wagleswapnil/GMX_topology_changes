import sys, os, re
import argparse
from protein_top_parser import parse_protein_top
from write_topology import write_protein_topology

def update_idx(i, idx, increment):
    if int(i) >= idx:
        i = str(int(i) + increment)
    return i

def shift_index(protein, idx, increment):
    idx = int(idx)
    increment = int(increment)
    if "atoms" in protein:
        for i, atom in enumerate(protein["atoms"]):
            if atom[0] == ";":
                continue
            elif atom[0] == "#":
                continue
            else:
                #nr, rest = atom.split(None, 1)
                #protein["atoms"][i] = update_idx(nr, idx, increment).rjust(6) + "    " + rest
                nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
                protein["atoms"][i] = update_idx(nr, idx, increment).rjust(6) + "  " + atype.rjust(6) + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "  " + aname.rjust(6) + "  " + update_idx(cgnr, idx, increment).rjust(6) + "  " + q.rjust(12) + "  " + m.rjust(10) + " " + rest

    if "bonds" in protein:
        for i, bond in enumerate(protein["bonds"]):
            if bond[0] == ";":
                continue
            else:
                ai, aj, rest = bond.split(None, 2)
                protein["bonds"][i] = update_idx(ai, idx, increment).rjust(6) + "    " + update_idx(aj, idx, increment).rjust(6) + "    " + rest
    if "pairs" in protein:
        for i, pair in enumerate(protein["pairs"]):
            if pair[0] == ";":
                continue
            else:
                ai, aj, rest = pair.split(None, 2)
                protein["pairs"][i] = update_idx(ai, idx, increment).rjust(6) + "    " + update_idx(aj, idx, increment).rjust(6) + "    " + rest
    if "angles" in protein:
        for i, angle in enumerate(protein["angles"]):
            if angle[0] == ";":
                continue
            else:
                ai, aj, ak, rest = angle.split(None, 3)
                protein["angles"][i] = update_idx(ai, idx, increment).rjust(6) + "    " + update_idx(aj, idx, increment).rjust(6) +"    " + update_idx(ak, idx, increment).rjust(6) + "    " + rest
    if "dihedrals" in protein:
        for i, dihedral in enumerate(protein["dihedrals"]):
            if dihedral[0] == ";":
                continue
            if dihedral[0] == "#":
                print("Comment in the dihedrals section. Stopping reindexing here.")
                break
            else:
                ai, aj, ak, al, rest = dihedral.split(None, 4)
                protein["dihedrals"][i] = update_idx(ai, idx, increment).rjust(6) + "    " + update_idx(aj, idx, increment).rjust(6) +"    " + update_idx(ak, idx, increment).rjust(6) + "    " + update_idx(al, idx, increment).rjust(6) + "    " + rest
    return protein

def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--topA", type=str, help="path for topology A, i.e., the topology of the wild type protein", required=True)
    parser.add_argument("--topB", type=str, help="path for topology B, i.e., the topology of the mutant protein", required=False)
    parser.add_argument("--idx", type=str, help="index, after which the addition/deletion is to be made", required=True)
    parser.add_argument("--increment", type=str, help="the amplitude of the addition/deletion", required=True)
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypes, protein = parse_protein_top(args.topA, "system1")
    protein = shift_index(protein, args.idx, args.increment)
    write_protein_topology(atomtypes, protein)
    return

if __name__== "__main__":
    main()
    
