import argparse
from read_topology import TopologyFile

def input_data():
    parser = argparse.ArgumentParser(description="Change GROMACS topologies")
    parser.add_argument("--top_path1", type=str, help="Path of the first TOP file")
    args = parser.parse_args()
    return args


def delete_residues(args):
    top_file1 = TopologyFile(args.top_path1)
    top_file1.read_topology()
    return


def main():
    args = input_data()
    delete_residues(args)
    return

if __name__ == "__main__":
    main()
