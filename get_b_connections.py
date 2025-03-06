import sys, os, re
import argparse
from protein_top_parser import parse_protein_top


def get_topB_indices(proteinB, resB):
    residues = {str(int(resB) - 1): "nminus1", resB: "n", str(int(resB) + 1): "nplus1"}
    bb_atoms = ["N", "CA", "C", "O", "H", "HA", "CB"]
    topB_indices = {}
    resB_indices = []
    resnameB = None
    for atom in proteinB["atoms"]:
        if atom[0] == ";" or atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if resnr == resB:
                topB_indices[nr] = residues[resnr] + aname
                resB_indices.append(nr)
                resnameB = residue
            elif resnr in residues and aname in bb_atoms:
                topB_indices[nr] = residues[resnr] + aname
    return topB_indices, resB_indices, resnameB

def get_resB_connections(proteinB, topB_bb_indices, resB_indices):
    bonds = []
    for i, bond in enumerate(proteinB["bonds"]):
        if bond[0] == ";" or bond[0] == "#":
            continue
        else:
            ai, aj, rest = bond.split(None, 2)
            if ai in resB_indices or aj in resB_indices:
                bonds.append(topB_bb_indices[ai] + "    " + topB_bb_indices[aj] + "    " + rest)
    proteinB["bonds"] = bonds

    pairs = []
    for i, pair in enumerate(proteinB["pairs"]):
        if pair[0] == ";" or pair[0] == "#":
            continue
        else:
            ai, aj, rest = pair.split(None, 2)
            if ai in resB_indices or aj in resB_indices:
                pairs.append(topB_bb_indices[ai] + "    " + topB_bb_indices[aj] + "    " + rest)
    proteinB["pairs"] = pairs

    angles = []
    for i, angle in enumerate(proteinB["angles"]):
        if angle[0] == ";" or angle[0] == "#":
            continue
        else:
            ai, aj, ak, rest = angle.split(None, 3)
            if ai in resB_indices or aj in resB_indices or ak in resB_indices:
                angles.append(topB_bb_indices[ai] + "    " + topB_bb_indices[aj] + "    " + topB_bb_indices[ak] + "    " + rest)
    proteinB["angles"] = angles

    dihedrals = []
    for i, dihedral in enumerate(proteinB["dihedrals"]):
        if dihedral[0] == ";" or dihedral[0] == "#":
            continue
        else:
            ai, aj, ak, al, rest = dihedral.split(None, 4)
            if ai in resB_indices or aj in resB_indices or ak in resB_indices or al in resB_indices:
                dihedrals.append(topB_bb_indices[ai] + "    " + topB_bb_indices[aj] + "    " + topB_bb_indices[ak] + "    " + topB_bb_indices[al] + "    " + rest)
    proteinB["dihedrals"] = dihedrals
    return proteinB

