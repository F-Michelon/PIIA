% INPUTS:
%   - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
%   - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
% OUTPUTS:  
%   - selgene(G): gene G is a selected gene belonging to the k genes
%   - affinity(I,J,X): a set of 4 cells where their expression for the k selected genes 
%     are the same for all cells. I, J, X and L belong to four different classes.

% Number of genes for a combination
#const k = 3.

% Generate combinations of k genes
{selgene(G) : pert(C,G,S,CL)} = k.

% Generate the corresponding perturbation vectors
selpert(E,V,S,C) :- selgene(V), pert(E,V,S,C).

% Generate an input predicate for genes that are inputs
selinput(G) :- selgene(G), input(G).

% Generate a equal/4 predicate for cells I and J, from two different classes, where their expression for the gene G is equal
equal(I,J,X,G) :- selpert(I,G,S1,C1), selpert(J,G,S2,C2), selpert(X,G,S3,C3), C1<C2, C2<C3, S1==S2, S2==S3.

% Count the number of equal/4 predicates between 3 cells of different classes, among the k selected genes
countequal(I,J,X,M) :- M={equal(I,J,X,_)}, selpert(I,_,_,C1), selpert(J,_,_,C2), selpert(X,_,_,C3), C1<C2, C2<C3.

% Generate, or not, an affinity/3 predicate when the number of countequal/4
% is equal to k for cells I and J of different classes, which guarantees that all genes are similarly expressed in both cells I and J
0{affinity(I,J,X)}1 :- countequal(I,J,X,k), I!=J, I!=X, J!=X, selpert(I,_,_,C1), selpert(J,_,_,C2), selpert(X,_,_,C3), C1<C2, C2<C3.

% Count the number of input genes expressed at 1 for each affinity/3
nbInputOnes(C,N) :- N={pert(C,G,1,_) : selinput(G)}, affinity(C,_,_).

% Forbid the affinity/3 where the number of input genes expressed at 1 is less than 1
:- affinity(C,_,_), nbInputOnes(C,N), N < 1.

% Generate a diff/3 predicate for cells I1 and I2, from a same class, where their expression for the gene G is different
diff(I1,I2,G) :- selpert(I1,G,S1,C1), selpert(I2,G,S2,C2), C1==C2, S1!=S2, I1<I2.

% Count the number of diff/3 predicates between 2 cells of the same class, among the k selected genes
countdiff(I1,I2,M) :- M={diff(I1,I2,_)}, I1<I2, selpert(I1,_,_,C1), selpert(I2,_,_,C2),C1==C2.

% Forbids countdiff/3 predicates for cells I1 and I2, from the first class
% where there is no difference in expression values for the k selected genes
:- countdiff(I1,I2,0), selpert(I1,_,_,C1), selpert(I2,_,_,C2), C1==C2, affinity(I1,_,_), affinity(I2,_,_), I1<I2.

% Make the same for the second class
:- countdiff(I1,I2,0), selpert(I1,_,_,C1), selpert(I2,_,_,C2), C1==C2, affinity(_,I1,_), affinity(_,I2,_), I1<I2.

% Make the same for the third class
:- countdiff(I1,I2,0), selpert(I1,_,_,C1), selpert(I2,_,_,C2), C1==C2, affinity(_,_,I1), affinity(_,_,I2), I1<I2.

% Maximize the number of pairs of cells
#maximize{1, I : affinity(I,_,_)}.

% Show selgene/1 and affinity/4 predicates
#show selgene/1.
#show affinity/3.