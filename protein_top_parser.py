import sys, os, re

def get_next_chunk(lines):
    count = 0
    section = None
    for line in lines:
        if line.strip().startswith("[") and line.strip().endswith("]"):
            break
        else:
            count += 1
    
    if line[0].startswith("[") and line[0].endswith("]"):
        section = line[1:-2].strip()
    chunk = lines[:count]
    lines[:count] = []

    return section, chunk, lines


def parse_protein_top(path, molname):
    topology={}
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    for i in range(0, 5):
        section, chunk, lines = get_next_chunk(lines)
        breakpoint()
        print(section, "\n".join(chunk))
    return



