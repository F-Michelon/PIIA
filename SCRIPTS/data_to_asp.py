import pandas as pd
import numpy as np

def read_file(pkn_nodes_file:str) ->list :
    with open(pkn_nodes_file, 'r') as f:
        return f.read().splitlines()

def load_data(filename):
    return pd.read_csv(filename, compression='zip')

B_path = '../matrices/toy_datasets/B.csv'
B = pd.read_csv(B_path)

pkn_nodes_file = '../pkn_data/nodes.txt'
pkn_no_predecessors_file = '../pkn_data/no_predecessors.txt'
pkn_no_successors_file = '../pkn_data/no_successors.txt'
genes_list = read_file(pkn_nodes_file)
inputs_list = read_file(pkn_no_predecessors_file)
readouts_list = read_file(pkn_no_successors_file)
inputs = []
readouts = []
intermediate = []
parents = []
for i, x in enumerate(inputs_list):
    if x in B.keys()[9:]:
        inputs.append(B.keys()[B.keys() == x][0])

for i, x in enumerate(readouts_list):
    if x in B.keys()[9:]:
        readouts.append(B.keys()[B.keys() == x][0])

for x in B.keys()[9:]:
    if x not in readouts and x not in inputs:
        intermediate.append(x)

for x in genes_list:
    tree = x.replace('/', ' ').split()
    for j in range(len(tree) - 1):
        if tree[j] in B.keys()[9:] and tree[j+1] in B.keys()[9:]:
            parents.append([B.keys()[B.keys() == tree[j]][0], B.keys()[B.keys() == tree[j+1]][0]])

"""
INPUTS:
 - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
 - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
"""

with open('../matrices/toy_datasets/asp_data.lp', 'w') as asp_data:
    for gene in inputs:
        asp_data.write(f"input({str(gene).lower().replace(' ', '_').replace('.', '')}).\n")
    for gene in readouts:
        asp_data.write(f"readout({str(gene).lower().replace(' ', '_').replace('.', '')}).\n")
    for gene in intermediate:
        asp_data.write(f"intermediate({str(gene).lower().replace(' ', '_').replace('.', '')}).\n")
    for p in parents:
        asp_data.write(f"parent({str(p[0]).lower().replace(' ', '_').replace('.', '')}, {str(p[1]).lower().replace(' ', '_').replace('.', '')}).\n")
    for index in B.index:
        for j in range(9,38):
            asp_data.write(f"pert({str(B.iloc[index,0]).lower().replace(' ', '_').replace('.', '')},{str(B.keys()[j]).lower().replace(' ', '_').replace('.', '')},{str(B.iloc[index,j]).lower().replace(' ', '_').replace('.', '')},{str(B.iloc[index,2]).lower().replace(' ', '_').replace('.', '')}).\n")
    asp_data.write("""

% Generate combinations of k genes
{selinput(G) : pert(C,G,S,CL), not intermediate(G)} = 1 .
{selinter(G): intermediate(G)} = 3-1.

1{selinput(G)  : parent(G,I), input(G)} :- selinter(I).

% Generate the corresponding perturbation vectors
selpert(E,V,S,C) :- selinput(V), pert(E,V,S,C).
selpert(E,V,S,C) :- selinter(V), pert(E,V,S,C).

% Generate a equal/3 predicate for cells I and J, from two different classes, where their expression for the gene G is equal
equal(I,J,G) :- selpert(I,G,S1,C1), selpert(J,G,S2,C2), I!=J, S1==S2.


pot_affi(I,J) :- k={equal(I,J,_)},  selpert(I,_,_,C1), selpert(J,_,_,C2), C1<C2, I!=J.

% Generate, or not, an affinity/2 predicate when the number of countequal/3
% is equal to k for cells I and J of different classes, which guarantees that all genes are similarly expressed in both cells I and J
0{affinity(I,J)}1 :- pot_affi(I,J).

% % Count the number of input genes expressed at 1 for each affinity/2
nbInputOnes(C,N) :- N={pert(C,G,1,_) : selinput(G), input(G)}, affinity(C,_).

% % Forbid the affinity/2 where the number of input genes expressed at 1 is less than 1
:- affinity(C,_), nbInputOnes(C,N), N < 1.

:- affinity(I,J1), affinity(I,J2), J1!=J2.
:- affinity(I1,J), affinity(I2,J), I1!=I2.

:- pot_affi(I1,J), pot_affi(I2,J), affinity(I2,J), I1<I2.
:- pot_affi(I,J1), pot_affi(I,J2), affinity(I,J2), J1<J2.


selgene(G) :- selinput(G).
selgene(G) :- selinter(G).

% Maximize the number of pairs of cells
#maximize{1, I : affinity(I,_)}.

% Show selgene/1 and affinity/2 predicates
#show selgene/1.
#show affinity/2.""")