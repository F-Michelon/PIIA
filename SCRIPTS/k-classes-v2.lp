% Search of the maximal set of N-uplets of cells where the expression 
% of the k genes selected are the same for the N cells.
%
% INPUTS:
%   - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
%   - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
% OUTPUTS:  
%   - selgene(G): gene G is a selected gene belonging to the k genes
%   - affinity(I,J,X): a pair of cells where their expression for the k selected genes 
%     are the same for both cells. I and J belong to two different classes.

% constante k
#const k = 3.

% les différentes classes et une référence
class(CL) :- pert(_,_,_,CL).
nb_class(N) :- N = #count {CL : class(CL)}.
ref(early_te).

% Generate combinations of k genes
{selgene(G) : pert(C,G,S,CL)} = k.

% la cellule C1 de la classe référence R a la même expression d'un gène parmi les gènes sélectionnés avec la cellule C2 d'une autre classe S
match_gene(C1,C2,CL,G) :- ref(R), pert(C1,G,E,R), pert(C2,G,E,CL), selgene(G), CL != R.
match_pert(C1,C2) :- match_gene(C1,C2,_,_). % ce prédicat simplifie les ligne suivantes

% les cellules C1 et C2 ont le même vecteur ie les k mêmes expressions de gènes
match_class(C1,C2,CL) :- k == #count {G : match_gene(C1,C2,CL,G)}, match_pert(C1,C2), class(CL).

% Un match est possible si on trouve le vecteur dans chaque classe
possible_match(C1) :- nb_class(N), N-1 == #count {CL : match_class(C1,_,CL)}, match_pert(C1,_).

% On associe ou non un match
0{match(C)}1 :- possible_match(C).

% C1 et C2 sont 2 cellules différentes dans la classe référence early_te
% S'il existe une cellule C3 dans un autre stage de développement qui a une affinité avec les 2 alors elles ont toutes le même vecteur
% C1<C2 permet de garder le plus petit C1 puis on supprime C2

meme_vecteur(C1, C2) :- C1 < C2, match_class(C1,C3,_), match_class(C2,C3,_). 
:- match(C2), meme_vecteur(C1,C2). 

% maximise les affinités
#maximize {1, C : match(C)}.

% renvoie les k gènes sélectionnés et le nombre d'affinité
#show selgene/1.
#show match/1.