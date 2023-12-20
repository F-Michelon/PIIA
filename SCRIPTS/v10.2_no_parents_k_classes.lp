% Search of the maximal set of pairs of cells where the expression 
% of the k genes selected are the same for the two cells of the pair.
%
% INPUTS:
%   - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
%   - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
% OUTPUTS:  
%   - selgene(G): gene G is a selected gene belonging to the k genes
%   - affinity(I,J): a pair of cells where their expression for the k selected genes 
%     are the same for both cells. I and J belong to two different classes.

#const k = 3.
#const i = 1.

% The x classes and the reference class
class(CL) :- pert(_,_,_,CL).
nb_class(N) :- N = #count {CL : class(CL)}.
ref(early_te).

% Generate combinations of k genes
{selgene(G) : pert(C,G,S,CL)} = k.

% Generate a equal/3 predicate for cells I and J, from two different classes, one being the reference classe R and the other a different classe S where their expression for the gene G is equal
equal(I,J,G,S) :- pert(I,G,S1,R), pert(J,G,S2,S), I!=J, S1==S2, S!=R.
match_pert(I,J) :- match_gene(I,J,_,_). %simplify the predicate


affi_pert(I,J,S) :- k == #count{G : equal(I,J,G,S)},  match_pert(I,J), class(S).

pot_affi(I) :- nb_class(N), N-1 == #count {S : affi_pert(I,_,S)}, match_pert(I,_).

% Generate, or not, an affinity/1 predicate when the number of countequal/3
% is equal to k for cells I and J of different classes, which guarantees that all genes are similarly expressed in both cells I and J
0{affinity(I)}1 :- pot_affi(I).

% I1 and I2 are 2 different cells in the reference class early_te
% If there is a J cell in another development stage which has an affinity with the 2, then they all have the same vector
:- affinity(I), affi_pert(I1,J,_), affi_perts(I2,J,_), I1 < I2.

% Maximize the number of pairs of cells
#maximize{1, I : affinity(I)}.

% Show selgene/1 and affinity/2 predicates
#show selgene/1.
#show affinity/2.
