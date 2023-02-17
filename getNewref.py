'''
Get new reference genome (root sequence) for a subtree.
'''
## Some code chunks taken from jmcbroome/pango-sequences/extract_strains.py 
# https://github.com/jmcbroome/pango-sequences/blob/0f72d87e3c07afd21d0d5b00ebb550f8d1c4dc8a/extract_strains.py#L26

# How to run:
'''
python3 getNewref.py \
  -t /home/shloka/data/public-latest.all.masked.pb.gz \
  -ref /home/shloka/data/NC_045512v2.fa \
  -c A.1 \
  -o /home/shloka/data/
'''

import bte

def parse_reference(refpath):
    refstr = []
    with open(refpath) as inf:
        for entry in inf:
            if entry[0] != ">":
                refstr.append(entry.strip())
    return "".join(refstr)

def process_mutstr(mstr):
    loc = int(mstr[1:-1])
    alt = mstr[-1]
    return (loc,alt)

def impute_haplotype(refstr, mutd):
    update = list(refstr)
    for m in mutd:
        loc,alt = process_mutstr(m)
        update[loc] = alt
    return "".join(update)


def main():
    refgenome = parse_reference(refpath = args.reference)
    newref = impute_haplotype(refgenome, mutd)
    of = open(f'{args.outputdir}newref{args.clade}.fa',"w+")
    print(f'>+{args.clade}',file=of)
    print(newref,file=of)
    of.close()
  
parser = argparse.ArgumentParser(
                      prog = 'getNewReference',
                      description = 'Get root sequence of the clade provided by user'
                      )
parser.add_argument('-t', '--tree', type=str, metavar='', required=True, help='SARS-CoV2 global phylogenetic tree; Protobuf file')
parser.add_argument('-ref', '--reference', type=str, metavar='', required=True, help='SARS-CoV2 reference genome sequence')
parser.add_argument('-c', '--clade', type=str, metavar='', required=True, help='clade name')
parser.add_argument('-o', '--outputdir', type=str, metavar='', required=True, help='new reference sequence path')
args = parser.parse_args()

mat = bte.MATree(args.tree)
clade_tree = mat.get_clade(args.clade)
root_node = clade_tree.root
mutd = clade_tree.get_haplotype(root_node.id)

main()
