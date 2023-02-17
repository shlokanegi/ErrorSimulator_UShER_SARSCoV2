from ete3 import Tree
import argparse

parser = argparse.ArgumentParser(
                    prog = 'Calculate RFD',
                    description = 'Calculate Robinson-Fould Distance between phastSim simulated "True Tree" and errorSimulator generated "Tree with errors"',
                    epilog = 'Provide trees in newick format'
                    )
parser.add_argument('-t1', '--tree1', type=str, metavar='', required=True, help='phastSim simulated tree/subtree in newick format')
parser.add_argument('-t2', '--tree2', type=str, metavar='', required=True, help='errorSimulator\'s tree with errors in newick format')

args = parser.parse_args()


with open(args.tree1, 'rb') as nwk1:
    tree1 = nwk1.read()
with open(args.tree1, 'rb') as nwk2:
    tree2 = nwk2.read()

t1 = Tree(args.tree1) # true tree
t2 = Tree(args.tree2) # simulated tree
rf = t1.robinson_foulds(t2)[0]
max_rf = t1.robinson_foulds(t2)[1]
common_leaves = t1.robinson_foulds(t2)[2]
parts_t1 = t1.robinson_foulds(t2)[3]
parts_t2 = t1.robinson_foulds(t2)[4] 
print(t1, t2)
print("RF distance is %s over a total of %s" %(rf, max_rf))
print("Partitions in tree2 that were not found in tree1:", parts_t1 - parts_t2)
print("Partitions in tree1 that were not found in tree2:", parts_t2 - parts_t1)

