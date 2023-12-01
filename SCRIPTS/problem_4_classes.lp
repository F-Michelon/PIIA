% Search of the maximal set of pairs of cells where the expression 
% of the k genes selected are the same for the two cells of the pair.
%
% INPUTS:
%   - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
%   - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
% OUTPUTS:  
%   - selgene(G): gene G is a selected gene belonging to the k genes
%   - affinity(I,J,K,L): a pair of cells where their expression for the k selected genes 
%     are the same for both cells. I and J belong to two different classes.

#const x = 3.

%constant for the number of classes
#const c = 4.

% Generate combinations of x genes
{selgene(G) : pert(C,G,S,CL)} = x .

% Generate the corresponding perturbation vectors
selpert(E,V,S,C) :- selgene(V), pert(E,V,S,C).

% Generate an input predicate for genes that are inputs
selinput(G) :- selgene(G), input(G).

% Generate a equal/5 predicate for cells I, J, K and L, from two different classes, where their expression for the gene G is equal
equal(I, J, K, L, G) :- selpert(I,G,S1,C1), selpert(J,G,S2,C2), selpert(K,G,S3,C3), selpert(L,G,S4,C4), C1<C2, C2<C3, C3<C4, S1==S2, S2==S3, S3==S4.

% Count the number of equal/5 predicates between 4 cells of different classes, among the x selected genes
countequal(I,J,K,L,M) :- M={equal(I,J,K,L,_)}, selpert(I,_,_,C1), selpert(J,_,_,C2), selpert(K,_,_,C3), selpert(L,_,_,C4), C1<C2, C2<C3, C3<C4.

% Generate, or not, an affinity/4 predicate when the number of countequal/5
% is equal to x for cells I, J, K and L of different classes, which guarantees that all genes are similarly expressed in both cells I and J
0{affinity(I, J, K, L)}1 :- countequal(I, J, K, L, x), I!=J, I!=K, I!=L, J!=K, J!=L, K!=L, selpert(I,_,_,C1), selpert(J,_,_,C2), selpert(K,_,_,C3), selpert(L,_,_,C4), C1<C2, C2<C3, C3<C4.

% Count the number of input genes expressed at 1 for each affinity/4
nbInputOnes(C,N) :- N={pert(C,G,1,_) : selinput(G)}, affinity(C,_,_,_).

% Forbid the affinity/4 where the number of input genes expressed at 1 is less than 1
:- affinity(C,_,_,_), nbInputOnes(C,N), N < 1.

% Generate a diff/3 predicate for cells I1 and I2, from a same class, where their expression for the gene G is different
diff(I1,I2,G) :- selpert(I1,G,S1,C1), selpert(I2,G,S2,C2), C1==C2, S1!=S2, I1<I2.

% Count the number of diff/3 predicates between 2 cells of the same class, among the x selected genes
countdiff(I1,I2,M) :- M={diff(I1,I2,_)}, I1<I2, selpert(I1,_,_,C1), selpert(I2,_,_,C2),C1==C2.

% Forbids countdiff/3 predicates for cells I1 and I2, from the first class
% where there is no difference in expression values for the x selected genes
:- countdiff(I1,I2,0), selpert(I1,_,_,C1), selpert(I2,_,_,C2), C1==C2, affinity(I1,_,_,_), affinity(I2,_,_,_), I1<I2.

% Make the same for the second class
:- countdiff(I1,I2,0), selpert(I1,_,_,C1), selpert(I2,_,_,C2), C1==C2, affinity(_,I1,_,_), affinity(_,I2,_,_), I1<I2.

% Make the same for the third class
:- countdiff(I1,I2,0), selpert(I1,_,_,C1), selpert(I2,_,_,C2), C1==C2, affinity(_,_,I1,_), affinity(_,_,I2,_), I1<I2.

% Make the same for the fourth class
:- countdiff(I1,I2,0), selpert(I1,_,_,C1), selpert(I2,_,_,C2), C1==C2, affinity(_,_,_,I1), affinity(_,_,_,I2), I1<I2.

% Maximize the number of pairs of cells
#maximize{1, I : affinity(I,_,_,_)}.

% Show selgene/1 and affinity/4 predicates
#show selgene/1.
#show affinity/4.
