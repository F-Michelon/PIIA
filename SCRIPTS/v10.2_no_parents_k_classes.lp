% Search of the maximal set of pairs of cells where the expression 
% of the k genes selected are the same for the two cells of the pair.
%
% INPUTS:
%   - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
% OUTPUTS:  
%   - selgene(G): gene G is a selected gene belonging to the k genes
%   - affinity(I): a pair of cells where their expression for the k selected genes 
%     are the same for both cells. I and J belong to two different classes.

#const k = 3.

% The N classes and the reference class
class(CL) :- pert(_,_,_,CL).
nb_class(N) :- N = #count {CL : class(CL)}.
ref(early_te).

% Generate combinations of k genes
{selgene(G) : pert(C,G,S,CL)} = k.

% Generate a equal/4 predicate for cells I and J, from two different classes, one being the reference class R and the other a different 
% class S where their expression for the gene G is equal
equal(I,J,G,S) :- ref(R), pert(I,G,S1,R), pert(J,G,S2,S), selgene(G), I!=J, S1==S2, S!=R.
match_pert(I,J) :- equal(I,J,_,_). %simplify the predicate

% Generate a affi_pert/3 predicate for cells I and J that have the same k expressions of genes G, ie they have the same boolean vector
affi_pert(I,J,S) :- k == #count{G : equal(I,J,G,S)},  match_pert(I,J), class(S).

% Generate a pot_affi/1 predicate where all the classes have the boolean vector
pot_affi(I) :- nb_class(N), N-1 == #count {S : affi_pert(I,_,S)}, match_pert(I,_).

% Generate, or not, an affinity/1 predicate where cell I of the reference class, has the same boolean vector than N-1 other cells from
% the other N-1 classes
0{affinity(I)}1 :- pot_affi(I).

% I1 and I2 are 2 different cells in the reference class
% If there is a cell J in another class which has an affinity with the 2, then they all have the same vector
:- affinity(I2), affi_pert(I1,J,_), affi_pert(I2,J,_), I1 < I2.

%comment
:- affinity(I), affi_pert(I,J1,S), affi_pert(I,J2,S), J1 < J2.

% Maximize the number of pairs of cells
#maximize{1, I : affinity(I)}.

% Show selgene/1 and affinity/1 predicates
#show selgene/1.
#show affinity/1.
