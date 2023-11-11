% Les atomes sont de la forme : 
% atom(numéro cellule, stage, nom gène, expression gène)

% Fichier des atomes
#include "données/B_3_classes.lp".

% constante k
#const k = 3.

% les différentes classes et une référence
class(S):- atom(_,S,_,_).
nb_class(N) :- N = #count {S : class(S)}.
ref(early_te).

% selectionne les k gènes 
k{selec_gene(G) : atom(_,_,G,_)}k.

% la cellule C1 de la classe référence R a la même expression d'un gène parmi les gènes sélectionnés avec la cellule C2 d'une autre classe S
match_gene(C1,C2,S,G) :- ref(R), atom(C1,R,G,E), atom(C2,S,G,E), selec_gene(G), S != R.
match_atom(C1,C2) :- match_gene(C1,C2,_,_). % ce prédicat simplifie les ligne suivantes
% les cellules C1 et C2 ont le même vecteur ie les k même expressions de gènes
match_stage(C1,C2,S) :- k == #count {G : match_gene(C1,C2,S,G)}, match_atom(C1,C2), class(S).
% Un match est possible si on trouve le vecteur dans chaque stage
possible_match(C1) :- nb_class(N), N-1 == #count {S : match_stage(C1,_,S)}, match_atom(C1,_).
% On associe ou non un match
0{match(C)}1 :- possible_match(C).

% C1 et C2 sont 2 cellules différentes dans la classe référence early_te
% S'il existe une cellule C3 dans un autre stage de développement qui a une affinité avec les 2 alors elles ont toutes le même vecteur
% C1<C2 permet de garder le plus petit C1 puis on supprime C2
meme_vecteur(C1, C2) :- C1 < C2, match_stage(C1,C3,_), match_stage(C2,C3,_). 
:- match(C2), meme_vecteur(C1,C2). 

% maximise les affinités
#maximize {1, C : match(C)}.

% renvoie les k gènes sélectionnés et le nombre d'affinité
#show selec_gene/1.
#show match/1.