import pandas as pd
import numpy as np
import random

def read_file(pkn_nodes_file:str) ->list :
    with open(pkn_nodes_file, 'r') as f:
        return f.read().splitlines()

def load_data(filename):
    return pd.read_csv(filename, compression='zip')

print("Write the file name you want to work on")
file_name = input()
path = '../DONNEES/toy_datasets/' + file_name + '.csv'
data = pd.read_csv(path)

pkn_no_predecessors_file = '../DONNEES/no_predecessors.txt'
pkn_no_successors_file = '../DONNEES/no_successors.txt'
inputs_list = read_file(pkn_no_predecessors_file)
readouts_list = read_file(pkn_no_successors_file)
inputs = []
readouts = []
intermediate = []
number_of_genes = len(data.keys()) - 9
for i, x in enumerate(inputs_list):
    if x in data.keys()[9:]:
        inputs.append(data.keys()[data.keys() == x][0])

for i, x in enumerate(readouts_list):
    if x in data.keys()[9:]:
        readouts.append(data.keys()[data.keys() == x][0])

for x in data.keys()[9:]:
    if x not in readouts and x not in inputs:
        intermediate.append(x)

# on crée artificiellement une moitié de gènes parents et fils pour le programme v10.2
permutation = random.sample(range(number_of_genes), number_of_genes)
sons_index = permutation[:len(permutation)//2]
parents_index = permutation[len(permutation)//2:]
parents = []
for i in range(len(sons_index)):
    parents.append([data.keys()[parents_index[i] + 9], data.keys()[sons_index[i] + 9]])

"""
INPUTS:
 - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
 - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
"""

with open('../DONNEES/toy_datasets/asp_data_' + file_name + '.lp', 'w') as asp_data:
    for gene in inputs:
        asp_data.write(f"input({str(gene).lower().replace(' ', '_').replace('.', '')}).\n")
    for gene in readouts:
        asp_data.write(f"readout({str(gene).lower().replace(' ', '_').replace('.', '')}).\n")
    for gene in intermediate:
        asp_data.write(f"intermediate({str(gene).lower().replace(' ', '_').replace('.', '')}).\n")
    for p in parents:
        asp_data.write(f"parent({str(p[0]).lower().replace(' ', '_').replace('.', '')}, {str(p[1]).lower().replace(' ', '_').replace('.', '')}).\n")
    for index in data.index:
        for j in range(9,38):
            asp_data.write(f"pert({str(data.iloc[index,0]).lower().replace(' ', '_').replace('.', '')},{str(data.keys()[j]).lower().replace(' ', '_').replace('.', '')},{str(data.iloc[index,j]).lower().replace(' ', '_').replace('.', '')},{str(data.iloc[index,2]).lower().replace(' ', '_').replace('.', '')}).\n")