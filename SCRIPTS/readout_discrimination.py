# Ce script développe des methodes de discrimination de readout pour du multi classe

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

class Discrimination:
    def __init__(self, data, score, optim, long_vect_bool):
        self.data = data
        self.score = score
        self.optimization = optim
        self.same_vect_bool = None
        self.bool_vect = None
        self.init = None
        self.long_vect_bool = long_vect_bool
        self.score = 0

    def find_same_vect_bool(self):
        """
        Fonction qui trouve les lignes du
        tableau pandas qui ont le même vecteur booléen
        
        Renvoie
        -------
        List_index_same_vect_bool : List[List[int]]
            Une liste contenant des sous listes.
            Chaque sous-liste correspond à un vecteur booléen différent
            Ces sous-listes contiennent les indexes (pour le dataframe)
            des cellules qui ont le même vecteur booléen
        List_bool_vect : List[List[bool]]
            Une liste contenant des sous listes.
            Chaque sous liste correspond à un vecteur booléen différent
        """
        N = len(self.data)
        List_bool_vect = []
        List_index_same_vect_bool = []
        for i in range(N):
            vect_bool = [[self.data.iloc[i][k] for k in range(self.long_vect_bool)]] #Le vecteur booléen de la ligne i du tableau
            if not(vect_bool in List_bool_vect): #S'il existe pas déja
                List_bool_vect.append(vect_bool) #On l'ajoute à la liste
                List_index_same_vect_bool.append([i]) #Et on ajoute son indice à la liste des numéros de chasue vecteurs boléens
            
            else : #S'il est dans la liste
                for j in range(len(List_bool_vect)): #On cherche quel vecteur c'est
                    if List_bool_vect[j] == vect_bool: #Si c'est le vecteur n° j de la liste
                        List_index_same_vect_bool[j].append(i) #On ajoute son indice à la liste
        
        self.same_vect_bool = List_index_same_vect_bool
        self.bool_vect = List_bool_vect

    def define_3_cel_for_each_vect(self):
        """
        Fonction qui, pour chaque vecteur booléen différent,
        prend au hasard trois cellules qui partagent ce même vecteur booléen
        mais qui appartiennent à des classes différentes
        
        Renvoie
        -------
        List_3_cels_for_each_vect : List[int]
            La liste des indexes dans le dataframe des cellules choisies,
            il y en a 12 au total, 3 cellules pour chaque vecteurs
            et il y a 4 vecteurs différents
        """
        List_3_cels_for_each_vect = []
        nb_vectors_diff = len(self.bool_vect)
        for v in range(nb_vectors_diff):
            
            List_3_cels = [0, 0, 0] #La liste des 3 cellules qu'on va choisir
            List_class_of_List_3_cel = [] #La liste des classes correspondantes à ces trois cellules
            add = 0
            
            while add != 3:
                #On prend une cellule aléatoirement dans la liste des cellules qui ont le même vecteur booléen 
                r = np.random.randint(0,len(self.same_vect_bool[v])) #On la choisi parmis la liste qu'on connait
                indice_cel = self.same_vect_bool[v][r] #On prend son indice dans le dataframe
                cel_class = self.data['classes'][indice_cel] #on prend aussi sa classe
                
                if not(cel_class in List_class_of_List_3_cel): #Si on a pas déja une cellule de cette classe
                    add += 1
                    List_3_cels[cel_class] = indice_cel #On l'ajoute à la liste
                    List_class_of_List_3_cel.append(cel_class) #On ajoute aussi sa classe à la liste
            
            #Une fois qu'on a choisi nos 3 cellules (une pour chaque classe)
            #On peut l'ajouter à la liste globale
            List_3_cels_for_each_vect.append(List_3_cels)
        
        self.init = List_3_cels_for_each_vect

    def choix3cel(self, index_vect_bool, cel_to_change=None):
        """
        Fonction qui prend en argument 3 cellules ayant le même vecteur booléen
        mais qui sont chacune de classe différentes.
        La fonction renvoie un triplet de cellules dont 2 sont les mêmes que celles
        prises en argument et la troisième fait partie de la même classe que celle
        qui a été supprimé du triplet précédent.
        """
        if cel_to_change is None:
            cel_to_change = np.random.randint(0, 3)
        index_cel_to_change = self.init[cel_to_change]
        class_to_change = self.data['classes'][index_cel_to_change]
        
        change = False
        while not(change):
            index_new_cel = random.sample(self.same_vect_bool[index_vect_bool], 1)[0]
            class_new_cel = self.data['classes'][index_new_cel]
            change = index_new_cel != index_cel_to_change and class_to_change == class_new_cel
        
        self.init[cel_to_change] = index_new_cel

    def calculate_score(self, gene, vecClass1, vecClass2):
        """
        Fonction qui, pour un gène,
        calcule le score de différence entre les vecteurs booléens.
        """
        return np.linalg.norm(self.data.iloc[vecClass1][gene].values - self.data.iloc[vecClass2][gene].values)

    def compute_score(self, gene):
        score = 0
        for i in range(3):
            for j in range(i + 1, 4):
                score += self.calculate_score(self, gene, self.init[i], self.init[j])
        return score / 6

    def compute_score_global(self):
        score = 0
        for gene in self.data.keys()[self.long_vect_bool:-2]:
            score += self.compute_score(gene)
        self.score = score / len(self.data.keys()[self.long_vect_bool:-2])

    def maximize_score(self):
        """
        Fonction qui cherche le maximum du score
        """

        self.find_same_vect_bool(self)
        self.define_3_cel_for_each_vect(self)
        iter = 0
        self.compute_score_global(self)
        while iter < 10:
            iter += 1
            init_copy = init.copy()
            for i in range(4):
                init_copy = init.copy()
                Discrimination.choix3cel(self, i, vect_find[0])
                new_score = Discrimination.compute_score_global(self, init_copy)
                if new_score > score:
                    init = init_copy.copy()
                    score = new_score
            if new_score > score:
                init = init_copy.copy()
            if iter%100 == 0: print(score, iter)
        return init

PATH = "../DONNEES/toy_datasets/readout_fictifs_D_3.csv"

D = pd.read_csv(PATH)

# #Test
discrimination = Discrimination(D, None, None, 120)
discrimination.maximize_score()
# choix3cel(0,test_init[0],test_find[0],D)
# print(test_init[0])
# print(test_init[1])
# print(calculate_score('139', test_init[0], test_init[1], D))
# print(compute_score('139', test_init, D))
# print(compute_score_global(test_init, D))


# import timeit
# temps_execution  = timeit.timeit("maximize_score(D)",globals=globals(),number = 1)
# print(temps_execution)

# courbe_recherche = [[0.745570745611005, 200],
#                     [0.7498832203326453, 300],
#                     [0.7504543307124744, 400],
#                     [0.7504543307124744, 500],
#                     [0.7504543307124744, 600],
#                     [0.7518329620768774, 700],
#                     [0.7518329620768774, 800],
#                     [0.7518329620768774, 900],
#                     [0.7518329620768774, 1000],
#                     [0.7518329620768774, 1100],
#                     [0.7518329620768774, 1200],
#                     [0.7518329620768774, 1300],
#                     [0.7518329620768774, 1400],
#                     [0.7518329620768774, 1500],
#                     [0.7518329620768774, 1600],
#                     [0.7540423627869316, 1700],
#                     [0.7540423627869316, 1800],
#                     [0.7586052478056338, 1900],
#                     [0.7586052478056338, 2000],
#                     [0.7586052478056338, 2100],
#                     [0.7586052478056338, 2200],
#                     [0.7586052478056338, 2300],
#                     [0.7586052478056338, 2400],
#                     [0.7586052478056338, 2500],
#                     [0.7586052478056338, 2600],
#                     [0.7586052478056338, 2700],
#                     [0.7586052478056338, 2800],
#                     [0.7586052478056338, 2900],
#                     [0.7586052478056338, 3000],
#                     [0.7586052478056338, 3100],
#                     [0.7586052478056338, 3200],
#                     [0.7586052478056338, 3300],
#                     [0.7695919195382717, 3400],
#                     [0.7695919195382717, 3500],
#                     [0.7695919195382717, 3600],
#                     [0.7695919195382717, 3700],
#                     [0.7695919195382717, 3800],
#                     [0.7695919195382717, 3900],
#                     [0.7695919195382717, 4000],
#                     [0.7695919195382717, 4100],
#                     [0.7695919195382717, 4200],
#                     [0.7695919195382717, 4300],
#                     [0.7695919195382717, 4400],
#                     [0.7695919195382717, 4500],
#                     [0.7695919195382717, 4600],
#                     [0.7724791976066105, 4700],
#                     [0.7724791976066105, 4800],
#                     [0.7724791976066105, 4900],
#                     [0.7724791976066105, 5000],
#                     [0.7724791976066105, 5100],
#                     [0.7724791976066105, 5200],
#                     [0.7724791976066105, 5300],
#                     [0.7724791976066105, 5400],
#                     [0.7724791976066105, 5500],
#                     [0.7724791976066105, 5600],
#                     [0.7724791976066105, 5700],
#                     [0.7724791976066105, 5800],
#                     [0.7724791976066105, 5900],
#                     [0.7724791976066105, 6000],
#                     [0.7724791976066105, 6100],
#                     [0.7724791976066105, 6200],
#                     [0.7724791976066105, 6300],
#                     [0.7724791976066105, 6400],
#                     [0.7724791976066105, 6500],
#                     [0.7724791976066105, 6600],
#                     [0.7724791976066105, 6700],
#                     [0.7724791976066105, 6800],
#                     [0.7724791976066105, 6900],
#                     [0.7724791976066105, 7000],
#                     [0.7724791976066105, 7100],
#                     [0.7724791976066105, 7200],
#                     [0.7724791976066105, 7300],
#                     [0.7724791976066105, 7400],
#                     [0.7724791976066105, 7500],
#                     [0.7724791976066105, 7600],
#                     [0.7724791976066105, 7700],
#                     [0.7724791976066105, 7800],
#                     [0.7724791976066105, 7900],
#                     [0.7724791976066105, 8000],
#                     [0.7724791976066105, 8100],
#                     [0.7724791976066105, 8200],
#                     [0.7724791976066105, 8300],
#                     [0.7724791976066105, 8400],
#                     [0.7724791976066105, 8500],
#                     [0.7724791976066105, 8600],
#                     [0.7724791976066105, 8700],
#                     [0.7724791976066105, 8800],
#                     [0.7724791976066105, 8900],
#                     [0.7724791976066105, 9000],
#                     [0.7724791976066105, 9100],
#                     [0.7724791976066105, 9200],
#                     [0.7724791976066105, 9300],
#                     [0.7724791976066105, 9400],
#                     [0.7724791976066105, 9500],
#                     [0.7724791976066105, 9600],
#                     [0.7724791976066105, 9700],
#                     [0.7724791976066105, 9800],
#                     [0.7724791976066105, 9900],
#                     [0.7724791976066105, 10000]]

# import matplotlib.pyplot as plt
# plt.figure()
# Y = [item[0] for item in courbe_recherche]
# X = [item[1] for item in courbe_recherche]
# plt.plot(X,Y)
# plt.ylabel("score")
# plt.xlabel("Nombre d'execution")
# plt.title("Évolution du score en fonction du nombre d'execution")
# plt.grid()
# plt.show()