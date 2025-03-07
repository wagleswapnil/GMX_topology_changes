#!/usr/bin/python
"""
Description: This script inserts sidechain atoms from a residue from topology B into a residue of topology A. It generates a B state (dummy) for the sidechain of residue of topology A, and an A state (dummy) for the sidechain of residue of topology B. For every atom added, the script updates the atomtypes list by adding the dummmy atomtype of the added atom. 

- The residue number of the sidechains of both topologies A (where the new atoms are inserted) and B (where the new residue atoms are taken from)) have to be the same for subsequent use (connecting these inserted atoms with precursor and subsequent amino acids) of the output topology. 

Written on: March 06th, 2025
Written by: Swapnil Wagle, swapnil.wagle92[at]gmail[.]com
"""

import sys, os, re
import argparse
from protein_top_parser import parse_protein_top
from write_topology import write_protein_topology

#Fuction to convert the list of atomtypes into a dictionary.
def get_atomtypes_dict(atomtypes):
    atomtypes_dict = {}
    for i, atomtype in enumerate(atomtypes):
        if atomtype[0] == ";":
            continue
        else:
            atomtypes_dict[atomtype.split()[0]] = atomtype.split(None,1)[1].strip()
    return atomtypes_dict


#Function to generate the B (dummy) state for the sidechain of residue from topology A. 
#The fuction also updates the atomtypes list with the newly added dummy atomtype.
def generate_Bstate(atomtypes, protein, A, backbone):
    atomtypes_dict = get_atomtypes_dict(atomtypes)
    bb_atoms = ["N", "C", "O", "CA", "H", "HA"]
    beginA = None
    endA = None
    if "-" in A:
        beginA = int(A.split("-")[0])
        endA = int(A.split("-")[1])
    else:
        beginA = int(A.strip())
        endA = int(A.strip())
    for i, atom in enumerate(protein["atoms"]):
        if atom[0] == ";":
            continue
        elif atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) >= beginA and int(nr) <= endA:
                if backbone == "no":
                    if aname in bb_atoms:
                        continue
                protein["atoms"][i] = nr.rjust(6) + "  " + atype.rjust(6) + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "  " + aname.rjust(6) + "  " + cgnr.rjust(6) + "  " + q.rjust(12) + "  " + m.rjust(10) + "    " + "DUM" + atype + "    " + "0.000000" + "    " + m.rjust(10) + " " + rest
                if "DUM" + atype in atomtypes_dict:
                    continue
                else:
                    atomtypes_dict["DUM" + atype] = "0.000000     0.000000   A     0.000000     0.000000"
                    atomtypes.append("DUM" + atype + "    " + "0.000000     0.000000   A     0.000000     0.000000\n")
    return atomtypes, protein


#Function to insert sidechain atoms for a residue from topology B. 
#The added atoms are interacting in state B, while state A for these each of these atoms is dummy. 
#A dummy atomtype for each of the added atom is also added to the atomtypes list. 
def generate_Astate(atomtypes, protein, B, backbone):
    atomtypes_dict = get_atomtypes_dict(atomtypes)
    bb_atoms = ["N", "C", "O", "CA", "H", "HA"]
    beginB = None
    endB = None
    if "-" in B:
        beginB = int(B.split("-")[0])
        endB = int(B.split("-")[1])
    else:
        beginB = int(B.strip())
        endB = int(B.strip())
    for i, atom in enumerate(protein["atoms"]):
        if atom[0] == ";":
            continue
        elif atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) >= beginB and int(nr) <= endB:
                if backbone == "no":
                    if aname in bb_atoms:
                        continue
                protein["atoms"][i] = nr.rjust(6) + "   " + "DUM" + atype + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "    " + "D" + aname + "  " + cgnr.rjust(6) + "  " + "0.000000" + "    " + m.rjust(10) + "    " + atype.rjust(6) + "    " + q.rjust(12) + "  " + m.rjust(10) + " " + rest
                if "DUM" + atype in atomtypes_dict:
                    continue
                else:
                    atomtypes_dict["DUM" + atype] = "0.000000     0.000000   A     0.000000     0.000000"
                    atomtypes.append("DUM" + atype + "    " + "0.000000     0.000000   A     0.000000     0.000000\n")
    return atomtypes, protein


#This function takes atoms from topology file B and inserts into the topology file A.
def merge_atoms(atomtypes, proteinA, idxA, proteinB, idxB, backbone):
    atomtypes_dict = get_atomtypes_dict(atomtypes)
    bb_atoms = ["N", "C", "O", "CA", "H", "HA"]
    bbB_charges = {}
    B_state_atoms = []
    beginA = None
    endA = None
    beginB = None
    endB = None
    if "-" in idxA:
        beginA = int(idxA.split("-")[0])
        endA = int(idxA.split("-")[1])
    else:
        beginA = int(idxA.strip())
        endA = int(idxA.strip())
    if "-" in idxB:
        beginB = int(idxB.split("-")[0])
        endB = int(idxB.split("-")[1])
    else:
        beginB = int(idxB.strip())
        endB = int(idxB.strip())

    insertion_idx = endA
    for i, atom in enumerate(proteinB["atoms"]):
        if atom[0] == ";":
            continue
        elif atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) >= beginB and int(nr) <= endB:
                if aname in bb_atoms:
                    bbB_charges[aname] = q
                else:
                    insertion_idx += 1
                    B_state_atoms.append(str(insertion_idx).rjust(6) + "   " + "DUM" + atype + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "    " + aname + "  " + str(insertion_idx).rjust(6) + "  " + "0.000000" + "    " + m.rjust(10) + "    " + atype.rjust(6) + "    " + q.rjust(12) + "  " + m.rjust(10) + " " + rest)
                    if "DUM" + atype in atomtypes_dict:
                        atomtypes_dict["DUM" + atype] = "0.000000     0.000000   A     0.000000     0.000000"
                        atomtypes.append("DUM" + atype + "    " + "0.000000     0.000000   A     0.000000     0.000000\n")
 

    insertion_point = 0
    for i, atom in enumerate(proteinA["atoms"]):
        if atom[0] == ";":
            continue
        elif atom[0] == "#":
            continue
        else:
            nr, atype, resnr, residue, aname, cgnr, q, m, rest = atom.split(None, 8)
            if int(nr) >= beginA and int(nr) <= endA:
                if aname in bb_atoms:
                    proteinA["atoms"][i] = nr.rjust(6) + "  " + atype.rjust(6) + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "  " + aname.rjust(6) + "  " + cgnr.rjust(6) + "  " + q.rjust(12) + "  " + m.rjust(10) + "    " + atype + "    " + bbB_charges[aname] + "    " + m.rjust(10) + " " + rest
                else:
                    proteinA["atoms"][i] = nr.rjust(6) + "  " + atype.rjust(6) + "  " + resnr.rjust(6) + "  " + residue.rjust(6) + "  " + aname.rjust(6) + "  " + cgnr.rjust(6) + "  " + q.rjust(12) + "  " + m.rjust(10) + "    " + "DUM" + atype + "    " + "0.000000" + "    " + m.rjust(10) + " " + rest
                    insertion_point = i
                if "DUM" + atype in atomtypes_dict:
                    continue
                else:
                    atomtypes_dict["DUM" + atype] = "0.000000     0.000000   A     0.000000     0.000000"
                    atomtypes.append("DUM" + atype + "    " + "0.000000     0.000000   A     0.000000     0.000000\n")

    proteinA["atoms"] = proteinA["atoms"][:insertion_point+1] + B_state_atoms + proteinA["atoms"][insertion_point+1:]
    return atomtypes, proteinA

#Function called by "check_atomtypes" to compare sigma and epsilon values
def compare_atomtype_values(value1, value2):
    sigma1, epsilon1 = value1.split()[-2:]
    sigma2, epsilon2 = value2.split()[-2:]
    if sigma1 ==  sigma2 and epsilon1 == epsilon2:
        return True
    else:
        return False


#Function to check if the atomtpes from topology files A and B are the same.
# For example, sigma and epsilon for an atomtype "C1" should be the same in the two files.
def check_atomtypes(atomtypesA, atomtypesB):
    atomtypes = []
    atomtypesA_dict = get_atomtypes_dict(atomtypesA)
    atomtypesB_dict = get_atomtypes_dict(atomtypesB)

    for key, value in atomtypesA_dict.items():
        if key in atomtypesB_dict:
            if compare_atomtype_values(value, atomtypesB_dict[key]):
                continue
            else:
                print("Error" + key + " is different between the two files")
                quit()
        else:
            print("Error" + "key" + " is absent from file2")

    for key, value in atomtypesB_dict.items():
        if key in atomtypesA_dict:
            if compare_atomtype_values(value, atomtypesA_dict[key]):
                continue
            else:
                print("Error" + key + " is different between the two files")
                quit()
        else:
            print("Error" + "key" + " is absent from file1")
    return

def input_data():
    parser = argparse.ArgumentParser(description="generate a dual topology file, where one residue of the protein is mutated into another")
    parser.add_argument("--topA", type=str, help="path for topology file A", required=True)
    parser.add_argument("--topB", type=str, help="path for topology file B", required=True)
    parser.add_argument("--A", type=str, help="range of indices, for which A state is interacting and B state is dummy", required=True)
    parser.add_argument("--B", type=str, help="range of indices, for which B state is interacting and A state is dummy", required=True)
    parser.add_argument("--backbone", type=str, help="should A and B states  be generated for the backbone, yes/no?", required=False, default="no")
    args = parser.parse_args()
    return args

def main():
    args = input_data()
    atomtypesA, proteinA = parse_protein_top(args.topA, "system1")
    atomtypesB, proteinB = parse_protein_top(args.topB, "system1")
    check_atomtypes(atomtypesA, atomtypesB)
    atomtypes, protein = merge_atoms(atomtypesA, proteinA, args.A, proteinB, args.B, args.backbone)
    write_protein_topology(atomtypes, protein)
    return

if __name__== "__main__":
    main()
    
