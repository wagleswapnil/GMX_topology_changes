import sys, os, re
import argparse
from protein_top_parser import parse_protein_top

def get_topA_indices(proteinA, resA, resnameB):
    residues = {str(int(resA) - 1): "nminus1", str(int(resA) + 1): "nplus1"}
    bb_atoms = ["N", "CA", "C", "O", "H", "HA", "CB"]
    real_bb_atoms = ["N", "CA", "C", "O", "H", "HA"]
    topA_indices = {}
    resA_indices = []
    nCB = None
    for atom in proteinA["atoms"]:
        if atom[0] == ";" or atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if resnr == resA and aname in real_bb_atoms:
                topA_indices["n" + aname] = nr
                resA_indices.append(nr)
            elif resnr == resA and residue == resnameB:
                topA_indices["n" + aname] = nr
                resA_indices.append(nr)
            elif resnr in residues and aname in bb_atoms:
                topA_indices[residues[resnr] + aname] = nr
                resA_indices.append(nr)
    return topA_indices, resA_indices

def make_resB_connections(proteinA, stateB_connections, topA_indices, resA_indices):
    print("".join(stateB_connections["bonds"]))
    for i, bond in enumerate(proteinA["bonds"]):
        if bond[0] == ";" or bond[0] == "#":
            continue
        else:
            ai, aj, rest = bond.split(None, 2)
            if ai in resA_indices and aj in resA_indices:
                for bondB in stateB_connections["bonds"]:
                    aiB, ajB, restB = bondB.split(None, 2)
                    if topA_indices[aiB] == ai and topA_indices[ajB] == aj:
                        proteinA["bonds"][i] = ai + "    " + aj + "    " + rest.split(';')[0].strip() + "  " + restB
                        stateB_connections["bonds"].remove(bondB)
                    elif topA_indices[aiB] == aj and topA_indices[ajB] == ai:
                        proteinA["bonds"][i] = ai + "    " + aj + "    " + rest.split(';')[0].strip() + "  "  + restB
                        stateB_connections["bonds"].remove(bondB)

    if len(stateB_connections["bonds"]) != 0:
        for i, bondB in enumerate(stateB_connections["bonds"]):
            aiB, ajB, restB = bondB.split(None, 2)
            if aiB in topA_indices and ajB in topA_indices:
                proteinA["bonds"].append(topA_indices[aiB] + "    " + topA_indices[ajB] + "    " + restB)


    for i, pair in enumerate(proteinA["pairs"]):
        if pair[0] == ";" or pair[0] == "#":
            continue
        else:
            ai, aj, rest = pair.split(None, 2)
            if ai in resA_indices and aj in resA_indices:
                for pairB in stateB_connections["pairs"]:
                    aiB, ajB, restB = pairB.split(None, 2)
                    if topA_indices[aiB] == ai and topA_indices[ajB] == aj:
                        proteinA["pairs"][i] = ai + "    " + aj + "    " + rest.split(';')[0].strip() + "  " + restB
                        stateB_connections["pairs"].remove(pairB)
                    elif topA_indices[aiB] == aj and topA_indices[ajB] == ai:
                        proteinA["pairs"][i] = ai + "    " + aj + "    " + rest.split(';')[0].strip() + "  " + restB
                        stateB_connections["pairs"].remove(pairB)
    
    if len(stateB_connections["pairs"]) != 0:
        for i, pairB in enumerate(stateB_connections["pairs"]):
            aiB, ajB, restB = pairB.split(None, 2)
            if aiB in topA_indices and ajB in topA_indices:
                proteinA["pairs"].append(topA_indices[aiB] + "    " + topA_indices[ajB] + "    " + restB)


    for i, angle in enumerate(proteinA["angles"]):
        if angle[0] == ";" or angle[0] == "#":
            continue
        else:
            ai, aj, ak, rest = angle.split(None, 3)
            if ai in resA_indices and aj in resA_indices and ak in resA_indices:
                for angleB in stateB_connections["angles"]:
                    aiB, ajB, akB, restB = angleB.split(None, 3)
                    if topA_indices[aiB] == ai and topA_indices[ajB] == aj and topA_indices[akB] == ak:
                        proteinA["angles"][i] = ai + "    " + aj + "    " + ak + "  " + rest.split(';')[0].strip() + "  " + restB
                        stateB_connections["angles"].remove(angleB)
                    elif topA_indices[aiB] == ak and topA_indices[ajB] == aj and topA_indices[akB] == ai:
                        proteinA["angles"][i] = ai + "    " + aj + "    " + ak + "  " + rest.split(';')[0].strip() + "  " + restB
                        stateB_connections["angles"].remove(angleB)

    if len(stateB_connections["angles"]) != 0:
        for i, angleB in enumerate(stateB_connections["angles"]):
            aiB, ajB, akB, restB = angleB.split(None, 3)
            proteinA["angles"].append(topA_indices[aiB] + "    " + topA_indices[ajB] + "    " + topA_indices[akB] + "    " + restB)


    for i, dihedral in enumerate(proteinA["dihedrals"]):
        if dihedral[0] == ";" or dihedral[0] == "#":
            continue
        else:
            ai, aj, ak, al, rest = dihedral.split(None, 4)
            if ai in resA_indices and aj in resA_indices and ak in resA_indices and al in resA_indices:
                for dihedralB in stateB_connections["dihedrals"]:
                    aiB, ajB, akB, alB, restB = dihedralB.split(None, 4)
                    if topA_indices[aiB] == ai and topA_indices[ajB] == aj and topA_indices[akB] == ak and topA_indices[alB] == al:
                        proteinA["dihedrals"][i] = ai + "    " + aj + "    " + ak + "    " + al + rest.split(';')[0].strip() + "  " + restB
                        stateB_connections["dihedrals"].remove(dihedralB)
                    elif topA_indices[aiB] == al and topA_indices[ajB] == ak and topA_indices[akB] == aj and topA_indices[alB] == ai:
                        proteinA["dihedrals"][i] = ai + "    " + aj + "    " + ak + "    " + al + rest.split(';')[0].strip() + "  " + restB
                        stateB_connections["dihedrals"].remove(dihedralB)

    if len(stateB_connections["dihedrals"]) != 0:
        for i, dihedralB in enumerate(stateB_connections["dihedrals"]):
            aiB, ajB, akB, alB, restB = dihedralB.split(None, 4)
            proteinA["dihedrals"].append(topA_indices[aiB] + "    " + topA_indices[ajB] + "    " + topA_indices[akB] + "    " + topA_indices[alB] + "    " + restB)
    return proteinA
