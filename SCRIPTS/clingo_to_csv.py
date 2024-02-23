import pandas as pd
import numpy as np
import warnings
import collections

warnings.filterwarnings("ignore")

def read_file(path:str) ->list :
    with open(path, 'r') as f:
        return f.read().splitlines()

output_path = '../RESULTATS/RESULTATS_NAUTILUS/K_CLASSES_WITHOUT_ZEROS/resultat_2023_PKN_earlyAndMediumAndLate_traite_2_v10.2_no_parents_k_classes-v2.txt'
output = read_file(output_path)

path = '../DONNEES/real_datasets/'
data = pd.read_csv(path + '2023_PKN_earlyAndMediumAndLate_traite_2.csv')

optimum = output[-2].split(' ')

setup = {"stimuli": [], "inhibitors": [], "readouts": ["PKM", "CEBPB", "ZEB1", "SOD1", "NKX3-1", "GSR", "HESX1", "DDIT3", "ETV5", "HEY2", "HEY1", "PSAT1", "CEBPD", "DEC1"]}
cells = []

for char in optimum:
    if 'selinput' in char:
        setup['stimuli'].append(data.columns[int(char.replace('selinput(','')[:-1]) + 8])
    elif 'selinter' in char:
        setup['inhibitors'].append(data.columns[int(char.replace('selinter(','')[:-1]) + 8])
    elif 'affinity' in char:
        cells.append(char.replace('affinity(','')[:-1])

# def table de hashage
for i, cell in enumerate(data['Name']):
    transformed_cell = str(cell).lower().replace(' ', '_').replace('.', '')
    if transformed_cell in cells:
        cells[cells.index(transformed_cell)] = cell

def to_classe(x):
    if x == 'early_TE':
        return 0
    elif x == 'medium_TE':
        return 1
    else:
        return 2

reduced_data = data[['Name', 'clusterUmap'] + setup['stimuli'] + setup['inhibitors']]

readouts = pd.read_csv(path + 'raw_mtx.csv', usecols=['Name'] + setup['readouts'])
readouts = readouts[readouts['Name'].isin(reduced_data['Name'])]
readouts = readouts.set_index('Name').reindex(reduced_data['Name'])
readouts['Name'] = readouts.index
readouts = readouts.set_index(reduced_data.index)
readouts = readouts[['Name'] + list(readouts.keys())[:-1]]

for readout in setup['readouts']:
    reduced_data[readout] = (2/np.pi) * np.arctan(readouts[readout])

# test
# for name in reduced_data['Name']:
#     for readout in setup['readouts']:
#         print(list(reduced_data[reduced_data['Name'] == name][readout]))
#         if list(reduced_data[reduced_data['Name'] == name][readout])[0] != list(readouts[readouts['Name'] == name][readout])[0]:
#             print(name)

reduced_data['classes'] = reduced_data['clusterUmap'].apply(to_classe)
reduced_data = reduced_data.drop(columns=['clusterUmap'])

print(reduced_data)

# on récupère tous les gènes ayant un vecteur booléens similaire aux cellules

# cell1 = reduced_data[reduced_data['Name'].isin(cells)]
# cell2 = reduced_data[reduced_data['Name'].isin(list(set(reduced_data['Name'].values) - set(cells)))]

cells_to_keep = []
vect_bools = reduced_data[reduced_data['Name'].isin(cells)][setup['stimuli'] + setup['inhibitors']].values

for index in reduced_data['Name']:
    vect_bool_index = reduced_data[reduced_data['Name'] == index][setup['stimuli'] + setup['inhibitors']].values[0]
    for vect_bool in vect_bools:
        add = True
        for i in range(len(vect_bool)):
            if vect_bool[i] != vect_bool_index[i]:
                add = False
        if add and index not in cells_to_keep:
            cells_to_keep.append(index)

reduced_data[reduced_data['Name'].isin(cells_to_keep)].to_csv(path + '2023_PKN_earlyAndMediumAndLate_traite_2_readouts.csv', index=False)