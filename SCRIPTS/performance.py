# Etude des performances
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join

path = "../RESULTATS/"

onlyfiles_with_zeros = [f for f in listdir(path + 'WITH_ZEROS') if isfile(join(path, f))]
onlyfiles_without_zeros = [f for f in listdir(path + 'WITHOUT_ZEROS') if isfile(join(path, f))]

program_used = []
dataset_used = []
success = []
answers = []
time = []
for file in onlyfiles_with_zeros:
    with open(path + file, 'r') as file:
        file_name = file.name.replace('/', ' ').split()[-1]
        info = file_name.replace('_', ' ').split()
        output = file.read().splitlines()
        
        program_used.append(info[-1][:-4])
        dataset_used.append(info[1] + ' with zeros')
        success.append('Answer' in output[3])
        if success[-1]:
            for i in range(len(output)):
                if 'OPTIMUM FOUND' in output[i]:
                    answers.append(output[i - 2].split())
        else:
            answers.append([])
        time.append(output[-2][15:20])

for file in onlyfiles_without_zeros:
    with open(path + file, 'r') as file:
        file_name = file.name.replace('/', ' ').split()[-1]
        info = file_name.replace('_', ' ').split()
        output = file.read().splitlines()
        
        program_used.append(info[-1][:-4])
        dataset_used.append(info[1]  + ' without zeros')
        success.append('Answer' in output[3])
        if success[-1]:
            for i in range(len(output)):
                if 'OPTIMUM FOUND' in output[i]:
                    answers.append(output[i - 2].split())
        else:
            answers.append([])
        time.append(output[-2][15:20])

performances = pd.DataFrame({'program_used': program_used,
                                'dataset_used': dataset_used, 
                                'success': success, 
                                'answers': answers, 
                                'time': time}).to_csv(path + 'performances')
    