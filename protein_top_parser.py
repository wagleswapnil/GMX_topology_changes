import sys, os, re

def remove_blank_lines(chunk):
    return [line for line in chunk if line.strip()]

def get_next_chunk(lines):
    count = 0
    section = None
    for line in lines:
        if line.strip().startswith("[") and line.strip().endswith("]") and count != 0:
            count += 1
            break
        elif line.strip().startswith("[") and line.strip().endswith("]") and count == 0:
            #section = line[1:-2].strip()
            section = line.replace('[', '').replace(']', '').strip()
            count += 1
        else:
            count += 1
    if section == None:
        section = "head_comment"
    chunk = remove_blank_lines(lines[1:count-1])
    lines[:count-1] = []
    return section, chunk, lines


#currently does not work if there are "bonds" or "angles" or "dihedrals" sections in the intermolecular interactions at the end of the file
#TODO: Only works if the topology file contains the topology of only the mutating moleculetype. Otherwise, it just keeps on reading the file and ultimately throws an error
def parse_protein_top(path, molname):
    #topology={ "head" : {}, "molecules": {}, "tail": {}}
    topology = {}
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    head_sections = ["head_comment", "defaults", "atomtypes", "nonbond_params", "other"]
    moleculetype_section = ["atoms", "bonds", "pairs", "angles", "dihedrals", "cmap", "virtual_sites2", "constraints", "exclusions", "other", "settles"]
    tail_sections = ["system", "molecules", "other"]
    main_section = None
    topology['head'] = {}
    topology['head']['atomtypes'] = []
    while len(lines) != 1:
        section, chunk, lines = get_next_chunk(lines)
        if section in head_sections:
            main_section = "head"
        elif section == "moleculetype":
            for l in chunk:
                if l.strip()[0] != ";":
                    main_section = l.split()[0]
                else:
                    continue
        elif section in tail_sections:
            main_section = "tail"
        topology.setdefault(main_section, {})
        topology[main_section][section] = chunk
    if molname in topology.keys():
        return topology['head']['atomtypes'], topology[molname]
    else:
        print(molname + " not in topology. Exiting!")
        quit()

