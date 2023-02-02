## Sampling Leaves using Error counts estimated by Poisson Distribution
import bte
import numpy as np
from scipy.stats import poisson

def parse_reference(refpath):
    refstr = []
    with open(refpath) as inf:
        for entry in inf:
            if entry[0] != ">":
                refstr.append(entry.strip())
    return "".join(refstr)


def dfs_traversal_for_error_addition(current_node, sample_leaf_nodes, previous_node_DNA_sequence, transition_matrix, base_to_idx_mapping, reference_genome):
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
    
    ERRORS
    11. current_ref_base          : reference allele which has to be changed
    12. current_alt_base          : alternate allele chosen randomly based on transition probability matrix
    13. transition_counts         : transition counts (later probability) corresponding the `current_ref_base`
  '''

  if(current_node.id == mat.root.id):
    ## current node is the root, so store the reference DNA sequence in `current_node_DNA_sequence`
    current_node_DNA_sequence = list(reference_genome)
  else:
    ## update `current_node_DNA_sequence` and `transition_matrix` here
      current_node_DNA_sequence = previous_node_DNA_sequence
      for mutation in current_node.mutations:
        mutation_index = int(mutation[1:-1])  # assuming 0-indexing
        ref_base = mutation[0]
        alt_base = mutation[-1]
        current_node_DNA_sequence[mutation_index] = alt_base     #### CHECK: Do we need to store sequences of internal nodes too?
        transition_matrix[base_to_idx_mapping[ref_base]][base_to_idx_mapping[alt_base]] += 1
  

  ## Incorporating errors into the sampled leaf nodes
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
            current_node_mutations.append(''.join([current_ref_base, str(error_site_idx), current_alt_base]))
            break

      print(f"Errors added: {current_node_mutations[-num_errors:]}")
      current_node.update_mutations(current_node_mutations)
      transition_matrix = transition_matrix + error_transition_matrix
      print(f"No. of mutations after error addition: {len(current_node.mutations)}")
      print()
        
  
  
  ## If current_node is not a leaf node
  else:
    for child_node in current_node.children:
      dfs_traversal_for_error_addition(child_node, sample_leaf_nodes, current_node_DNA_sequence, transition_matrix, base_to_idx_mapping, reference_genome)

  
  ## update `current_node_DNA_sequence` and `transition_matrix` here
  for mutation in current_node.mutations:
    mutation_index = int(mutation[1:-1])  # assuming 0-indexing
    ref_base = mutation[0]
    alt_base = mutation[-1]
    current_node_DNA_sequence[mutation_index] = ref_base
    transition_matrix[base_to_idx_mapping[ref_base]][base_to_idx_mapping[alt_base]] -= 1

  return


def main():
  reference_genome = parse_reference(refpath = '/home/shloka/data/newref.fa') # input from argparser
  node_list = mat.depth_first_expansion()
  tree_mutation_count = sum([len(node.mutations) for node in node_list])
  expected_error_count = (error_rate * tree_mutation_count)/100
  error_count = poisson.rvs(expected_error_count)

  print(f"Total number of errors being incorporated: {error_count}\n")
  leaf_node_list = mat.get_leaves()
  leaf_mutation_count_dict = {}

  for leaf_node in leaf_node_list:
    leaf_mutation_count_dict[leaf_node.id] = len(mat.get_haplotype(leaf_node.id))

  # Sampling the nodes using the probability distribution, and have created `sample_leaf_nodes`
  leaf_ids = np.array(list(leaf_mutation_count_dict.keys()))
  leaf_mutation_counts = np.array(list(leaf_mutation_count_dict.values()))
  sum_leaf_mutations = sum(leaf_mutation_counts)

  sample_leaf_ids = list(np.random.choice(leaf_ids, size=error_count))
  sample_leaf_nodes = {x:sample_leaf_ids.count(x) for x in sample_leaf_ids}
  print(f"Total leaves on tree: {len(leaf_node_list)}")
  print(f"No. of leaves with error addition: {len(sample_leaf_nodes)}\n")

  sequence = []
  base_to_idx_mapping = {'A': 0, 'T': 1, 'C': 2, 'G': 3 }
  transition_matrix = np.ones((4, 4))
  dfs_traversal_for_error_addition(mat.root, sample_leaf_nodes, sequence, transition_matrix, base_to_idx_mapping, reference_genome)  
  
  return 

mat = bte.MATree(pb_file = '/home/shloka/data/phastSim_output/sars-cov-2_simulation_output.mat.pb')
error_rate = 0.1 # in percentage (will be taken as input)
main()
