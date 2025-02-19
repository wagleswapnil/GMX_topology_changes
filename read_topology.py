import sys, os
class TopologyFile():
    def __init__(self, top_path):
        self.top_path = top_path
        self.topology = dict()

    def determine_chunk(self, lines):
        print(lines[0])
        return


    def read_topology(self):
        sections = ["moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals" ]
        with open(self.top_path, 'r') as f:
            for line in f:
                if "moleculetype" in line:
                    print(f.readline())
        f.close()
        return self.topology
    
def main():
    topology = TopologyFile(sys.argv[1])
    print(topology.read_topology())

if __name__ == "__main__":
    main()

