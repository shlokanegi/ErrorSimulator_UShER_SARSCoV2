'''
Author  : Shloka Negi, shnegi@ucsc.edu
Purpose : Error Simulator which allows the addition of random errors, reversions or amplicon dropouts with contamination on simulated real trees
Inputs  : MAT.pb (Mutation Annotated Tree Object)
        : Reference genome (or root sequence, if using a subtree)
        : random error rate - error count calculated through a poisson draw, with mu = total mutations on tree * random error rate
        : reversion error rate - error count calculated through a binomial draw, with p = reversion error rate and n = total mutations on tree
        : amplicon dropout rate - no. of amplicon dropouts calculated through a poisson draw, with mu = no. of amplicon primers in ARTIC4.1 * amplicon dropout rate
Outputs : VCF file with real mutations + errors
        : metadata.TSV file for visualizing erroneous trees on Taxonium (CoV2Tree.org)
Usage   : python3 errorSimulator.py -t tree.mat.pb -ref reference.fa -r 0.5 -rev 0.04 -ad 0.8

'''

import bte
import numpy as np
from scipy.stats import poisson
from scipy.stats import binom
import argparse

def parse_reference(refpath):
    refstr = []
    with open(refpath) as inf:
        for entry in inf:
            if entry[0] != ">":
                refstr.append(entry.strip())
    return [*refstr[0]]


def chromosome_update_to_mutations(node_list):
    '''
    Fixing the chromosome issue in the existing mutations of the tree
    '''

    for current_node in node_list:
        mutations = current_node.mutations #[A2G, T3462G, ...., etc.]
        mutation_with_chromosome_list = ['NC_045512v2:'+mutation for mutation in mutations]
        current_node.update_mutations(mutation_with_chromosome_list)
    return


def reconstruct_sequence_slice(leaf_node, reference_genome, amplicon_range, mutation_list):
    '''
    helper method for amplicon_dropout()
    '''
    sequence_slice = reference_genome[amplicon_range[0]: amplicon_range[1]]   
    for mutation in mutation_list:
        index = int(mutation[1:-1])
        alt_base = mutation[-1]
        if (index >= amplicon_range[0]) and (index < amplicon_range[1]):
            sequence_slice[index-amplicon_range[0]] = alt_base
    return sequence_slice


def dfs_traversal_and_error_addition(current_node, leaf_ids, sample_leaf_nodes, previous_node_DNA_sequence, transition_matrix, base_to_idx_mapping, reference_genome):
    '''
    Description: Follows DFS to traverse each node in the tree.

    Variables used:

    1. current_node               : current node the execution is at (function called from main using root node for the first time).
    2. current_node_DNA_sequence  : DNA sequence at the current_node.
    3. previous_node_DNA_sequence : DNA sequence of the previous node from which the current_node was generated.
    4. reference_genome           : Ideally sequence of reference genome (here - sequence at root) --> list of characters

    5. mutation_index             : location on the sequence where mutation occurred.
    6. ref_base                   : Reference allele
    7. alt_base                   : Alternate allele
    8. transition_matrix          : 2D numpy array with counts of all transitions and trasversion corresponding to the current_node.

    9. num_errors                 : number of errors to introduce in the sampled_leaf_node
    10. error_site_idx            : Possible error site chosen at random -> int

    RANDOM ERRORS
    11. current_ref_base          : reference allele which has to be changed
    12. current_alt_base          : alternate allele chosen randomly based on transition probability matrix
    13. transition_counts         : transition counts (later probability) corresponding the `current_ref_base`
    '''

    if(current_node.id == mat.root.id):
        ## current node is the root, so store the reference DNA sequence in `current_node_DNA_sequence`
        current_node_DNA_sequence = reference_genome
    else:
        ## update `current_node_DNA_sequence` and `transition_matrix` here
        current_node_DNA_sequence = previous_node_DNA_sequence
        for mutation in current_node.mutations:
            mutation_list = mutation.split(":")
            mutation = mutation_list[-1]

            mutation_index = int(mutation[1:-1])  # assuming 0-indexing
            ref_base = mutation[0]
            alt_base = mutation[-1]
            current_node_DNA_sequence[mutation_index] = alt_base     #### CHECK: Do we need to store sequences of internal nodes too?
            transition_matrix[base_to_idx_mapping[ref_base]][base_to_idx_mapping[alt_base]] += 1


    ## Incorporating random errors into the sampled leaf nodes
    if current_node.id in sample_leaf_nodes:
        num_errors = sample_leaf_nodes[current_node.id]
        if(num_errors != 0):
            print(f"No. of mutations for nodeID {current_node.id} is {len(current_node.mutations)}.\nNo. of errors to be incorporated: {num_errors}")
            current_node_mutations = current_node.mutations

            error_transition_matrix = np.zeros((4, 4))
            for i in range(num_errors): # Randomly incorporate `num_errors` errors in current_node
                while True:
                    error_site_idx = np.random.randint(len(reference_genome))   # choose error site `error_site_idx`
                    if current_node_DNA_sequence[error_site_idx] == reference_genome[error_site_idx]:   # means this site is non-mutated
                        current_ref_base = current_node_DNA_sequence[error_site_idx]
                        ## Choosing current_alt_base based on transition probability matrix
                        transition_counts = transition_matrix[base_to_idx_mapping[current_ref_base]]
                        transition_counts = transition_counts / transition_counts.sum()
                        current_alt_base = np.random.choice(['A', 'T', 'C', 'G'], p=transition_counts)   # we have sampled the `current_alt_base` using the transition probability matrix
                        error_transition_matrix[base_to_idx_mapping[current_ref_base]][base_to_idx_mapping[current_alt_base]] += 1
                        current_node_mutations.append(''.join(["NC_045512v2:", current_ref_base, str(error_site_idx), current_alt_base]))
                        # current_node_mutations.append(''.join([current_ref_base, str(error_site_idx), current_alt_base]))
                        break

            print(f"Errors added: {current_node_mutations[-num_errors:]}")
            # print("Non-updated mutation list: ", current_node.mutations)
            current_node.update_mutations(current_node_mutations)
            # print("Updated mutation list: ", current_node.mutations)
            transition_matrix = transition_matrix + error_transition_matrix
            print(f"No. of mutations after error addition: {len(current_node.mutations)}")
            print()
        
   
  ## If current_node is not a leaf node
    else:
        for child_node in current_node.children:
            dfs_traversal_and_error_addition(child_node, leaf_ids, sample_leaf_nodes, current_node_DNA_sequence, transition_matrix, base_to_idx_mapping, reference_genome)

    
    ## update `current_node_DNA_sequence` and `transition_matrix` here
    for mutation in current_node.mutations:
        mutation_list = mutation.split(":")
        mutation = mutation_list[-1]

        mutation_index = int(mutation[1:-1])  # assuming 0-indexing
        ref_base = mutation[0]
        alt_base = mutation[-1]
        current_node_DNA_sequence[mutation_index] = ref_base
        transition_matrix[base_to_idx_mapping[ref_base]][base_to_idx_mapping[alt_base]] -= 1
    return


def reversion_addition(leaf_node_list, reversion_count, leaf_mutation_count_dict):
    '''
    Adds reversions to all leaf nodes, with number of reversions per leaf node chosen through a binomial draw (at each leaf)
    '''

    leaf_mutation_probability_distribution = np.array([leaf_mutation_count_dict[leaf_node.id] for leaf_node in leaf_node_list])
    leaf_mutation_probability_distribution = leaf_mutation_probability_distribution / leaf_mutation_probability_distribution.sum()
    leaf_reversion_sample_list = list(np.random.choice(leaf_node_list, p=leaf_mutation_probability_distribution, size=reversion_count))
    sample_leaf_node_id_dict = {x.id:leaf_reversion_sample_list.count(x) for x in leaf_reversion_sample_list}
    samples_with_reversions_dict = {leaf : ("Yes" if leaf in sample_leaf_node_id_dict else "No") for leaf in leaf_mutation_count_dict}

    for leaf_node in leaf_node_list:
        if leaf_node.id in sample_leaf_node_id_dict:
            mutations_root_to_leaf = list(mat.get_haplotype(leaf_node.id))
            num_reversions_current_leaf = sample_leaf_node_id_dict[leaf_node.id]
            reversion_indices_list = list(np.random.randint(len(mutations_root_to_leaf), size=num_reversions_current_leaf))
            print(f"No. of reversions to be added on this leaf: {len(reversion_indices_list)}")
            current_leaf_mutations = leaf_node.mutations
            print(f"No. of mutations of leaf before adding reversions: {len(current_leaf_mutations)}")
            for idx in reversion_indices_list:
                mutation = mutations_root_to_leaf[idx]
                current_ref_base = mutation[0]
                error_site_idx = int(mutation[1:-1])
                current_alt_base = mutation[-1]
                current_leaf_mutations.append(''.join(["NC_045512v2:", current_alt_base, str(error_site_idx), current_ref_base]))
            
            print(f"No. of mutations of leaf after adding reversions: {len(current_leaf_mutations)}")
            print()
            leaf_node.update_mutations(current_leaf_mutations)
    
    return samples_with_reversions_dict


def amplicon_dropout(amplicon_ranges_list, n_amplicon, amplicon_dropout_count, leaf_node_list, reference_genome, leaf_mutation_count_dict):
    '''
    Using a real primer schene (ARTIC4.1) and an amplicon_dropout_rate (provided by the user), this function "drops out" an amplicon and 
    replaces it with variation from elsewhere on the tree 

    Source leaf      : Leaf where the amplicon is dropped out and replaced by the amplicon sequence from replacement leaf corresponding to same amplicon range 
    Replacement leaf : Leaf chosen for replacing the dropped out amplicon on source leaf with variation of itself.
    '''

    sampled_source_leaf_list = list(np.random.choice(leaf_node_list, size=amplicon_dropout_count))
    sampled_source_leaf_dict = {x:sampled_source_leaf_list.count(x) for x in sampled_source_leaf_list}
    sampled_amplicons_dict_list = {x:list(np.random.choice([i for i in range(n_amplicon)], size=sampled_source_leaf_dict[x], replace=False)) for x in sampled_source_leaf_dict}
    sampled_amplicons_dict_list = {key:[amplicon_ranges_list[x] for x in val] for key,val in sampled_amplicons_dict_list.items()}
    sampled_replacement_leaf_dict_list = {}
    for source_leaf, cnt in sampled_source_leaf_dict.items():
        rep_leaf_list = []
        for i in range(cnt):
            while True:
                rep_leaf = np.random.choice(leaf_node_list)
                if rep_leaf.id != source_leaf.id:
                    break
            rep_leaf_list.append(rep_leaf)
        sampled_replacement_leaf_dict_list[source_leaf] = rep_leaf_list

    # sampled_source_leaf_list = [A, B, C, A, B, A]
    # sampled_source_leaf_dict = {A: 3, B: 2, C: 1}
    # sampled_amplicons_dict_list = {A: [(50, 100), (80, 120), (170, 190)], B: [(70, 200), (140, 160)], C: [(210, 235)]}
    # sampled_replacement_leaf_dict_list = {A: [C, B, C], B: [C, A], C: [A]}

    samples_with_amplicon_drops_dict = {leaf : "No" for leaf in leaf_mutation_count_dict}

    for source_leaf in sampled_source_leaf_dict:
        for i in range(sampled_source_leaf_dict[source_leaf]):
            replacement_leaf = sampled_replacement_leaf_dict_list[source_leaf][i]
            amplicon_range = sampled_amplicons_dict_list[source_leaf][i]
            source_mutation_list = list(mat.get_haplotype(source_leaf.id))
            replacement_mutation_list = list(mat.get_haplotype(replacement_leaf.id))
            source_sequence_slice = reconstruct_sequence_slice(source_leaf, reference_genome, amplicon_range, source_mutation_list)
            replacement_sequence_slice = reconstruct_sequence_slice(replacement_leaf, reference_genome, amplicon_range, replacement_mutation_list)
            refgenome_sequence_slice = reference_genome[amplicon_range[0]: amplicon_range[1]]
            
            source_reversion_list = []
            for idx in range(amplicon_range[1]-amplicon_range[0]):
                if (source_sequence_slice[idx] == replacement_sequence_slice[idx]) and source_sequence_slice[idx] != refgenome_sequence_slice[idx] :
                    print("same mutation present in both source and replacement")
                if (source_sequence_slice[idx] != refgenome_sequence_slice[idx]) and (source_sequence_slice[idx] != replacement_sequence_slice[idx]):  # revert here!!
                    # first performing reversion in source_leaf (deleting mutations)
                    source_reversion_list.append(source_sequence_slice[idx] + str(idx+amplicon_range[0]) + refgenome_sequence_slice[idx])

            source_addition_list = []
            for idx in range(amplicon_range[1]-amplicon_range[0]):
                if replacement_sequence_slice[idx] != refgenome_sequence_slice[idx] and (source_sequence_slice[idx] != replacement_sequence_slice[idx]):  # add here!!
                    # second performing addition in source_leaf (from replacement_leaf)
                    source_addition_list.append(refgenome_sequence_slice[idx] + str(idx+amplicon_range[0]) + replacement_sequence_slice[idx])
            
            if len(source_reversion_list) != 0 or len(source_addition_list) != 0:
                samples_with_amplicon_drops_dict[source_leaf.id] = "Yes"
                print(f"SOURCE mutations = {source_mutation_list}")
                print(f"REPLACEMENT mutations = {replacement_mutation_list}")
                print(f"Number of mutations dropped from source leaf lying in amplicon range {amplicon_range} is {len(source_reversion_list)}")
                print(f"Number of mutations added to the source leaf from replacement leaf lying in amplicon range {amplicon_range} is {len(source_addition_list)}")
                print()

                # now add the mutations to the source and replacement leaf's mutation lists
                source_mutation_list.extend(source_reversion_list)
                source_mutation_list.extend(source_addition_list)
                
                # updating the mutations back into nodes
                source_leaf.update_mutations(source_mutation_list)
                print(f"Final number of mutations on source leaf: {len(source_leaf.mutations)}")
                print()

    return samples_with_amplicon_drops_dict



def main():
    reference_genome = parse_reference(refpath = args.reference)
    node_list = mat.depth_first_expansion()

    ## Updating tree with CHROM:mutation
    chromosome_update_to_mutations(node_list)

    tree_mutation_count = sum([len(node.mutations) for node in node_list])
    leaf_node_list = mat.get_leaves()
    print(f"Total leaves on tree: {len(leaf_node_list)}")

    leaf_mutation_count_dict = {}

    for leaf_node in leaf_node_list:
        leaf_mutation_count_dict[leaf_node.id] = len(mat.get_haplotype(leaf_node.id))
    
    ### Read all amplicon ranges from the BED file
    with open('/home/shloka/data/amplicon/SARS-CoV-2.insert.bed', 'r') as f:
      amplicon_ranges_list = [(int(line.split()[1]), int(line.split()[2])) for line in f.readlines()]
    n_amplicon = len(amplicon_ranges_list)

    reversion_count = binom.rvs(n=tree_mutation_count, p=reversion_error_rate)
    expected_error_count = (error_rate * tree_mutation_count)
    error_count = poisson.rvs(expected_error_count)
    expected_dropout_count = (amplicon_dropout_rate * n_amplicon)
    amplicon_dropout_count = poisson.rvs(expected_dropout_count)

    ## Initializing dictionaries for writing metadata file 
    samples_with_random_errors_dict = samples_with_reversions_dict = samples_with_amplicon_drops_dict = {leaf : "No" for leaf in leaf_mutation_count_dict}


    ########### REVERSIONS #############
    if reversion_count != 0:
        print(f"Total number of reversions to be added on the tree: {reversion_count}")
        samples_with_reversions_dict = reversion_addition(leaf_node_list, reversion_count, leaf_mutation_count_dict)

    
    ########### RANDOM ERRORS #############
    if error_count != 0:
        # Sampling the nodes for random error addition, and have created `sample_leaf_nodes`
        leaf_ids = np.array(list(leaf_mutation_count_dict.keys()))
        sample_leaf_ids = list(np.random.choice(leaf_ids, size=error_count))
        sample_leaf_nodes = {x:sample_leaf_ids.count(x) for x in sample_leaf_ids}
        
        
        samples_with_random_errors_dict = {leaf : ("Yes" if leaf in sample_leaf_nodes else "No") for leaf in leaf_mutation_count_dict}
        # {node1: "Yes", node2: "No", etc....}


        print(f"Total number of random errors to be incorporated on the tree: {error_count}")
        print(f"No. of leaves with random error addition: {len(sample_leaf_nodes)}\n")
        
        sequence = []
        base_to_idx_mapping = {'A': 0, 'T': 1, 'C': 2, 'G': 3 }
        transition_matrix = np.ones((4, 4))
        dfs_traversal_and_error_addition(mat.root, leaf_ids, sample_leaf_nodes, sequence, transition_matrix, base_to_idx_mapping, reference_genome)  

    ############ FIRST: Make sure amplicon_ranges_list has been defined ###########
    try:
        if amplicon_ranges_list == None:
            raise ValueError("")
    except ValueError:
        print("Define amplicon_ranges_list FIRST !!!")

    ## Run amplicon-dropout method
    if amplicon_dropout_rate != 0:
        print(f"{amplicon_dropout_count} amplicon dropouts will be made out of {n_amplicon} possible inserts.\n")
        samples_with_amplicon_drops_dict = amplicon_dropout(amplicon_ranges_list, n_amplicon, amplicon_dropout_count, leaf_node_list, reference_genome, leaf_mutation_count_dict)

    #Write into VCF - with original mutations + errors + reversions + amplicon dropouts
    mat.write_vcf(vcf_file = "subtree_errors.vcf") 

    # Write a metadata file for taxonium
    with open("/home/shloka/data/taxonium/metadata_errors.tsv", "w") as m:
        m.write("strain" + "\t" + "random_errors" + "\t" + "reversions" + "\t" + "amplicon_dropouts" + "\n")
        for leaf in samples_with_random_errors_dict:
            m.write(leaf + "\t" + samples_with_random_errors_dict[leaf] + "\t" + samples_with_reversions_dict[leaf] + "\t" + samples_with_amplicon_drops_dict[leaf] + "\n")

    return


parser = argparse.ArgumentParser(
                    prog = 'ErrorSimulator',
                    description = 'Simulate random errors, reversions and amplicon dropout in SARS-CoV2 phylogenetic tree',
                    epilog = 'Provide each error rate, else it would default to 0'
                    )
parser.add_argument('-t', '--tree', type=str, metavar='', required=True, help='phastSim output tree/subtree protobuf file')
parser.add_argument('-ref', '--reference', type=str, metavar='', required=True, help='root/reference genome sequence')
parser.add_argument('-r', '--random', type=float, metavar='', nargs='?', const=0, default=0, help='Random error rate, DEFAULT: 0')
parser.add_argument('-rev', '--reversion', type=float, metavar='', nargs='?', const=0, default=0, help='Reversion rate, DEFAULT: 0')
parser.add_argument('-ad', '--amplicon_drop', type=float, metavar='', nargs='?', const=0, default=0, help='Amplicon dropout rate, DEFAULT: 0')
args = parser.parse_args()

mat = bte.MATree(args.tree)
error_rate = args.random
reversion_error_rate = args.reversion
amplicon_dropout_rate = args.amplicon_drop
main()

