% Search of the maximal set of pairs of cells where the expression 
% of the k genes selected are the same for the two cells of the pair.
%
% INPUTS:
%   - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
%   - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
% OUTPUTS:  
%   - selgene(G): gene G is a selected gene belonging to the k genes
%   - affinity(I,J,X): a pair of cells where their expression for the k selected genes 
%     are the same for both cells. I and J belong to two different classes.

#const k = 3.
#const i = 1.
% Generate combinations of k genes
{selinput(G) : pert(C,G,S,CL), not intermediate(G)} = i .
{selinter(G): intermediate(G)} = k-i.

%1{selinput(G)  : parent(G,I), input(G)} :- selinter(I).

% Generate the corresponding perturbation vectors
selpert(E,V,S,C) :- selinput(V), pert(E,V,S,C).
selpert(E,V,S,C) :- selinter(V), pert(E,V,S,C).

% Generate a equal/4 predicate for cells I and J, from two different classes, where their expression for the gene G is equal
equal(I,J,X,G) :- selpert(I,G,S1,C1), selpert(J,G,S2,C2), selpert(X,G,S3,C3), I!=J, I!=X, J!=X, S1==S2, S2==S3.

pot_affi(I,J,X) :- k={equal(I,J,X,_)},  selpert(I,_,_,C1), selpert(J,_,_,C2), selpert(X,_,_,C3), C1<C2, C2<C3, I!=J, I!=X, J!=X.

% Generate, or not, an affinity/3 predicate when the number of countequal/4
% is equal to k for cells I and J of different classes, which guarantees that all genes are similarly expressed in both cells I and J
0{affinity(I,J,X)}1 :- pot_affi(I,J,X).

% % Count the number of input genes expressed at 1 for each affinity/3
nbInputOnes(C,N) :- N={pert(C,G,1,_) : selinput(G), input(G)}, affinity(C,_,_).

% % Forbid the affinity/3 where the number of input genes expressed at 1 is less than 1
:- affinity(C,_,_), nbInputOnes(C,N), N < 1.

:- affinity(I,J1,_), affinity(I,J2,_), J1!=J2.
:- affinity(I1,J,_), affinity(I2,J,_), I1!=I2.

:- affinity(_,J,X1), affinity(_,J,X2), X1!=X2.
:- affinity(_,J1,X), affinity(_,J2,X), J1!=J2.

:- pot_affi(I,J1,_), pot_affi(I,J2,_), affinity(I,J2,_), J1<J2.
:- pot_affi(I1,J,_), pot_affi(I2,J,_), affinity(I2,J,_), I1<I2.

:- pot_affi(_,J,X1), pot_affi(_,J,X2), affinity(_,J,X2), X1<X2.
:- pot_affi(_,J1,X), pot_affi(_,J2,X), affinity(_,J2,X), J1<J2.

selgene(G) :- selinput(G).
selgene(G) :- selinter(G).

% Maximize the number of pairs of cells
#maximize{1, I : affinity(I,_,_)}.

% Show selgene/1 and affinity/3 predicates
#show selgene/1.
#show affinity/3.