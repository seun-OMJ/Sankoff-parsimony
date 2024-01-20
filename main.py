from Bio import Phylo
import numpy as np
import sys
from Bio import SeqIO
import itertools

given_tree = sys.argv[1]
given_seq = sys.argv[2]

sequences = {record.id: str(record.seq) for record in SeqIO.parse(given_seq, "fasta")}

tree = Phylo.read(given_tree, 'newick')
bases = ['A', 'C', 'G', 'T']    

# function to calculate substitution costs
def substitution_cost(base1, base2):
    if base1 == base2:
        return 0
    elif (base1 in 'AG' and base2 in 'AG') or (base1 in 'CT' and base2 in 'CT'):
        return 1  # Transition
    else:
        return 2  # Transversion

root = tree.root
n = len(sequences[list(sequences.keys())[0]])  
cost_matrices = {}
for node in tree.find_clades(order = 'postorder'): 
    
    #initialize the leaf nodes matrices   
    if node.is_terminal():
        cost_matrices[node] = np.zeros((4, n))
        idx = list(sequences.keys()).index(node.name)
       
        for i, base in enumerate(sequences[node.name]):
            for j, b in enumerate(bases):
                cost = 0 if base == b else float('inf')
                cost_matrices[node][j][i] = cost
                
    # fill the internal nodes score matrices          
    else:
        cost_matrices[node] = np.zeros((4, n))
        for i in range(n):
            for j in range(4):
                min_cost0 = float('inf')
                min_cost1 = float('inf')
                for k in range(4):
                    temp_cost0 = (cost_matrices[node.clades[0]][k][i] + substitution_cost('ACGT'[j], 'ACGT'[k]))
                    temp_cost1 = (cost_matrices[node.clades[1]][k][i] + substitution_cost('ACGT'[j], 'ACGT'[k]))
                    if temp_cost0 < min_cost0:
                        min_cost0 = temp_cost0
                    if temp_cost1 < min_cost1:
                        min_cost1 = temp_cost1
                cost_matrices[node][j][i] = min_cost1 + min_cost0

# Output cost matrix and optimal sequences for the root
root_cost_matrix = cost_matrices[root]
optimal_sequences = []
for i in range(4):
    sequence = ''
    for j in range(n):
        sequence += 'ACGT'[np.argmin(root_cost_matrix[i])]
    optimal_sequences.append(sequence)


for i, base in enumerate(bases):
    print(f"{base}: {' '.join(map(str, root_cost_matrix[i]))}")

transposed_matrix = list(zip(*cost_matrices[root]))
min_costs = [min(column) for column in transposed_matrix]
min_indices = [[i for i, cost in enumerate(row) if cost == min_cost] for row, min_cost in zip(transposed_matrix, min_costs)]
strings = []
for indices in itertools.product(*min_indices):
    string = ''.join(bases[idx] for idx in indices)
    strings.append(string)

for string in strings:
    print(string)