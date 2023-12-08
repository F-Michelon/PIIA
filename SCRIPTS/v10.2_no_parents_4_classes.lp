% Search of the maximal set of pairs of cells where the expression 
% of the k genes selected are the same for the two cells of the pair.
%
% INPUTS:
%   - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
%   - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
% OUTPUTS:  
%   - selgene(G): gene G is a selected gene belonging to the x genes
%   - affinity(I,J,K,L): a pair of cells where their expression for the x selected genes 
%     are the same for both cells. I and J belong to two different classes.

#const x = 3.
#const i = 1.
% Generate combinations of x genes
{selinput(G) : pert(C,G,S,CL), not intermediate(G)} = i .
{selinter(G): intermediate(G)} = x-i.

%1{selinput(G)  : parent(G,I), input(G)} :- selinter(I).

% Generate the corresponding perturbation vectors
selpert(E,V,S,C) :- selinput(V), pert(E,V,S,C).
selpert(E,V,S,C) :- selinter(V), pert(E,V,S,C).

% Generate a equal/5 predicate for cells I and J, from two different classes, where their expression for the gene G is equal
equal(I, J, K, L, G) :- selpert(I,G,S1,C1), selpert(J,G,S2,C2), selpert(K,G,S3,C3), selpert(L,G,S4,C4), C1<C2, C2<C3, C3<C4, S1==S2, S2==S3, S3==S4.


pot_affi(I, J, K, L) :- x={equal(I,J,K,L,_)},  selpert(I,_,_,C1), selpert(J,_,_,C2), selpert(K,_,_,C3), selpert(L,_,_,C4), C1<C2, C2<C3, C3<C4, I!=J, I!=K, I!=L, J!=K, J!=L, K!=L.

% Generate, or not, an affinity/4 predicate when the number of countequal/5
% is equal to x for cells I and J of different classes, which guarantees that all genes are similarly expressed in both cells I and J
0{affinity(I,J, K, L)}1 :- pot_affi(I,J, K, L).

% % Count the number of input genes expressed at 1 for each affinity/4
nbInputOnes(C,N) :- N={pert(C,G,1,_) : selinput(G), input(G)}, affinity(C,_,_,_).

% % Forbid the affinity/4 where the number of input genes expressed at 1 is less than 1
:- affinity(C,_,_,_), nbInputOnes(C,N), N < 1.

:- affinity(I,J1,_,_), affinity(I,J2,_,_), J1!=J2.
:- affinity(I1,J,_,_), affinity(I2,J,_,_), I1!=I2.

:- affinity(_,J,K3,_), affinity(_,J,K4,_), K3!=K4.
:- affinity(_,J1,K,_), affinity(_,J2,K,_), J1!=J2.

:- affinity(_,_,K,L1), affinity(_,_,K,L2), L1!=L2.
:- affinity(_,_,K5,L), affinity(_,_,K6,L), K5!=K6.

:- pot_affi(I1,J,_,_), pot_affi(I2,J,_,_), affinity(I2,J,_,_), I1<I2.
:- pot_affi(I,J1,_,_), pot_affi(I,J2,_,_), affinity(I,J2,_,_), J1<J2.

:- pot_affi(_,J1,K,_), pot_affi(_,J2,K,_), affinity(_,J2,K,_), J1<J2.
:- pot_affi(_,J,K7,_), pot_affi(_,J,K8,_), affinity(_,J,K8,_), K7<K8.

:- pot_affi(_,_,K9,L), pot_affi(_,_,K10,L), affinity(_,_,K10,L), K9<K10.
:- pot_affi(_,_,K,L1), pot_affi(_,_,K,L2), affinity(_,_,K,L2), L1<L2.


selgene(G) :- selinput(G).
selgene(G) :- selinter(G).

% Maximize the number of pairs of cells
#maximize{1, I : affinity(I,_,_,_)}.

% Show selgene/1 and affinity/4 predicates
#show selgene/1.
#show affinity/4.
