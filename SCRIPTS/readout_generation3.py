'''
Ce script génère le fichier MIDAS nécessaire aux calculs de réseaux booléens sur caspo.
A chacun de l'ajouter dans son répertoire DONNEES mais sans le versionné sur git ! Donc fini les "git add --all". 
'''

import pandas as pd
import json
import numpy as np

PATH = "../DONNEES/real_datasets/raw_mtx.csv"

with open('../DONNEES/setup.json') as setup_file:
    setup = json.load(setup_file)
with open('../DONNEES/perturbations.json') as perturbation_file:
    perturbation = json.load(perturbation_file)

D = pd.read_csv(PATH, usecols=['Name', 'clusterUmap'] + setup["stimuli"] + setup["inhibitors"] + setup['readouts'])[['Name', 'clusterUmap'] + setup["stimuli"] + setup["inhibitors"] + setup['readouts']]
D = D.set_index('Name').loc[list(np.array(perturbation).flatten())].sort_values(by=['clusterUmap'], ascending=False)

time_cols = [f'DA:{name}' for name in setup['readouts']]
D[time_cols] = 10 * np.ones((len(D) ,len(time_cols)))
print(D)
D = D.drop(columns='clusterUmap')

D[D[setup["stimuli"] + setup["inhibitors"]] == 1] = 0
D[D[setup["stimuli"] + setup["inhibitors"]] > 1] = 1
D[setup["readouts"]] = (2/np.pi) * np.arctan(D[setup['readouts']])

dict_cols = {}
cols = []
for name in setup['stimuli']:
    dict_cols[name] = f'TR:{name}'
    cols.append(f'TR:{name}')
for name in setup['inhibitors']:
    dict_cols[name] = f'TR:{name}i'
    cols.append(f'TR:{name}i')
cols = cols + time_cols
readouts = []
for name in setup['readouts']:
    dict_cols[name] = f'DV:{name}'
    cols.append(f'DV:{name}')
    readouts.append(f'DV:{name}')

D = D.rename(columns=dict_cols)
D['TR:Toy:CellLine'] = 1
cols = ['TR:Toy:CellLine'] + cols
D = D[cols]

D_0 = D.iloc[:len(D)//2].copy()
D_0[readouts] = 0
D_0[time_cols] = 0
D_1 = D.iloc[len(D)//2:].copy()
D_1[readouts] = 0
D_1[time_cols] = 0
D_medium = pd.concat([D_0,D.iloc[:len(D)//2]], axis=0)
D_late = pd.concat([D_1,D.iloc[len(D)//2:]], axis=0)
D_late.to_csv("../DONNEES/dataset_late_midas .csv", index=False)
D_medium.to_csv("../DONNEES/dataset_medium_midas .csv", index=False)
print(D_0)
print(D_late)
print(D_medium)