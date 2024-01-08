# Ce script développe des methodes de discrimination de readout pour du multi classe

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

# TODO décrire la tête des données formattage
# paramétriser l'ensemble des données numérique

class Discrimination:
    def __init__(self, data, optim, long_vect_bool):
        self.data = data
        self.score = 0
        self.optimization = optim
        self.same_vect_bool = None
        self.bool_vect = None
        self.init = None
        self.long_vect_bool = long_vect_bool
        self.nb_classe = len(self.data['classes'].unique())

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
                index_bool_classes = [[] for i in range(self.nb_classe)]
                index_bool_classes[self.class_cell(i)].append(i)
                List_index_same_vect_bool.append(index_bool_classes) #Et on ajoute son indice à la liste des numéros de chaque vecteurs boléens
            
            else : #S'il est dans la liste
                for j in range(len(List_bool_vect)): #On cherche quel vecteur c'est
                    if List_bool_vect[j] == vect_bool: #Si c'est le vecteur n° j de la liste
                        List_index_same_vect_bool[j][self.class_cell(i)].append(i) #On ajoute son indice à la liste
        
        self.same_vect_bool = List_index_same_vect_bool
        self.bool_vect = List_bool_vect

    def class_cell(self, index_cell):
        # Renvoie le numéro de la classe de la cellule
        return self.data['classes'][index_cell]

    def define_3_cel_for_each_vect(self):
        """
        Fonction qui, pour chaque vecteur booléen différent,
        prend au hasard trois cellules qui partagent ce même vecteur booléen
        mais qui appartiennent à des classes différentes
        
        Renvoie
        -------
        List_3_cels_for_each_vect : List[List[int]]
            La liste des indexes dans le dataframe des cellules choisies,
            il y en a 12 au total, 3 cellules pour chaque vecteurs
            et il y a 4 vecteurs différents
        """
        List_3_cels_for_each_vect = []
        nb_vectors_diff = len(self.bool_vect)
        for v in range(nb_vectors_diff):
            
            List_3_cels = [0, 0, 0] #La liste des 3 cellules qu'on va choisir
            i = 0
            
            while i < self.nb_classe:
                #On prend une cellule aléatoirement dans la liste des cellules qui ont le même vecteur booléen 
                indice_cel = random.sample(self.same_vect_bool[v][i], 1)[0] #On prend son indice dans le dataframe
                List_3_cels[i] = indice_cel
                i = i + 1
                
            #Une fois qu'on a choisi nos 3 cellules (une pour chaque classe)
            #On peut l'ajouter à la liste globale
            List_3_cels_for_each_vect.append(List_3_cels)
        
        self.init = List_3_cels_for_each_vect

    def choix3cel(self, index_vect_bool, nb_cell_to_change=1):
        ### TODO ajout argument pour changer plus de 1 cellules
        
        """
        Fonction qui prend en argument 3 cellules ayant le même vecteur booléen
        mais qui sont chacune de classe différentes.
        La fonction renvoie un triplet de cellules dont 2 sont les mêmes que celles
        prises en argument et la troisième fait partie de la même classe que celle
        qui a été supprimé du triplet précédent.
        
        nb_cell_to_change = Nombre de cellules à changer, supérieur à 1, FIXME add max value depend nb cell per bool vect
        """
        
        cells_to_change = np.random.randint(0, self.nb_classe, nb_cell_to_change)
        for cell_to_change in cells_to_change:
            index_cel_to_change = self.init[index_vect_bool][cell_to_change]
            class_to_change = self.class_cell(index_cel_to_change)
            
            change = False
            while not(change):
                index_new_cel = random.sample(self.same_vect_bool[index_vect_bool][class_to_change], 1)[0]
                change = index_new_cel != index_cel_to_change
            
            self.init[index_vect_bool][cell_to_change] = index_new_cel

    def calculate_score(self, id_gene, vecClass1, vecClass2):
        """
        Fonction qui, pour un gène,
        calcule le score de différence entre les vecteurs booléens.
        """
        return np.linalg.norm(self.data.iloc[vecClass1][id_gene].values - self.data.iloc[vecClass2][id_gene].values)

    def compute_score(self, gene):
        score = 0
        for i in range(self.nb_classe):
            for j in range(i + 1, len(self.bool_vect)):
                score += self.calculate_score(gene, self.init[i], self.init[j])
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
        self.find_same_vect_bool()
        self.define_3_cel_for_each_vect()
        iter = 0
        no_change_iter = 0
        self.compute_score_global()
        while iter < 1000 and no_change_iter < 100:
            iter += 1
            old_init = self.init.copy()
            old_score = self.score.copy()
            for i in range(len(self.bool_vect)):
                self.choix3cel(i, nb_cell_to_change=int(1 + (no_change_iter > 10) + (no_change_iter > 30)))
                self.compute_score_global()
                if old_score >= self.score:
                    self.init = old_init
                    self.score = old_score
                    no_change = True
                else:
                    no_change_iter = 0
                    no_change = False
            no_change_iter += 1 * (no_change)
            if iter%10 == 0: print(self.score, iter, no_change_iter)

PATH = "../DONNEES/toy_datasets/readout_fictifs_D_3.csv"

D = pd.read_csv(PATH)

# #Test
discrimination = Discrimination(D, None, 120)
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