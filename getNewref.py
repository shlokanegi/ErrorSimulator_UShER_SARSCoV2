'''
Author  : Shloka Negi, shnegi@ucsc.edu; some code chunks taken from jmcbroome/pango-sequences/extract_strains.py 
Purpose : Get new reference genome (root sequence) for a subtree.
Inputs  : output VCFs from hap.py comparison
Usage   : python3 getNewref.py \
            -t public-latest.all.masked.pb.gz \
            -ref NC_045512v2.fa \
            -c B.1.1 \
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
