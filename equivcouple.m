% equivcouple   Combination of equivalent spins 
%
%  [F,N] = equivcouple(I,n)
%
%  The states due to n spins-I can be combined to give
%  a set of independent spins. Their quantum numbers
%  are returned in F, their respective abundances in N.
%
%  Example:
%
%     [F,N] = equivcouple(1/2,5)
%
%  5 spins-1/2 give rise to a first-order splitting pattern
%  [1 5 10 10 5 1] (see the function equivsplit). This can be
%  decomposed into one spin-5/2, four spin-3/2 and five spin-1/2
%  according to
%
%          5  5          5 spins-1/2
%       4  4  4  4       4 spins-3/2
%    1  1  1  1  1  1    1 spin-5/2
%   ------------------
%    1  5 10 10  5  1    sum
%
%  so F = [2.5 1.5 0.5] and N = [1 4 5].
%
%  In group theoretical terms, this corresponds to the reduction
%  of a product of n irreps of dimension 2*I of the rotation group
%  into a direct sum of irreps (Clebsch-Gordan decomposition).

% see J.H.Freed, G.K.Fraenkel, J.Chem.Phys. 39(2), 326-348 (1963)
% eq. (4.32)

function [F,N] = equivcouple(I,n)

if (nargin==0), help(mfilename); return; end

% Special case n=1: no action
if (n==1), F = I; N = n; return; end

% (1) List of reduced spin quantum numbers
F = I*n:-1:0; 

% (2) Compute number of spins for each reduced spin quantum number
E = ones(1,2*I+1);
RowVec = 1;
for q = 1:n
  RowVec = conv(RowVec,E);
end
RowVec = [0, RowVec];
dRowVec = diff(RowVec);

N = dRowVec(1:length(F));

% Mathematical basis:
% Reduction of tensor product representations of rotation group
% using Clebsch-Gordan direct sum decomposition. For two spins:
% D^(j1)xD(j2) = sum_{j=|j1-j2|}^{j1+j2} D^(j)
% For multiple spins, recursive.

return
