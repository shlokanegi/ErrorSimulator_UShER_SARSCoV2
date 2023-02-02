'''
Get new reference genome (root sequence) for a subtree.
'''

import bte

mat = bte.MATree(pb_file = '/home/shloka/data/public-latest.all.masked.pb.gz')
clade_tree = mat.get_clade('A.2.2')
root_node = clade_tree.root
mutd = clade_tree.get_haplotype(root_node.id)

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

refgenome = parse_reference('/home/shloka/data/NC_045512v2.fa')
newref = impute_haplotype(refgenome, mutd)
of = open('/home/shloka/data/newref.fa',"w+")
print(">"+'A.2.2',file=of)
print(newref,file=of)
of.close()
