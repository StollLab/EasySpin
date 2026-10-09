% equivcouple   Combination of equivalent spins
%
%  [F,N] = equivcouple(I,n)
%
%  The states due to n spins-I can be combined to give a set of independent
%  spins. Their quantum numbers are returned in F, their respective
%  abundances in N.
%
%  Input:
%    I    spin quantum number of the equivalent spins (0, 1/2, 1, 3/2, ...)
%    n    number of equivalent spins (0, 1, 2, ...)
%
%  Output:
%    F    row vector of spin quantum numbers of the coupled spins,
%         in decreasing order
%    N    row vector of abundances, one for each element of F
%
%  Example:
%
%     [F,N] = equivcouple(1/2,5)
%
%  5 spins-1/2 give rise to a first-order splitting pattern [1 5 10 10 5 1]
%  (see the function equivsplit). This can be decomposed into one spin-5/2,
%  four spin-3/2 and five spin-1/2 according to
%
%          5  5          5 spins-1/2
%       4  4  4  4       4 spins-3/2
%    1  1  1  1  1  1    1 spin-5/2
%   ------------------
%    1  5 10 10  5  1    sum
%
%  so F = [2.5 1.5 0.5] and N = [1 4 5].
%
%  In group theoretical terms, this corresponds to the reduction of a
%  product of n irreps of dimension 2*I+1 of the rotation group into a
%  direct sum of irreps (Clebsch-Gordan decomposition).
%
%  The results are exact as long as all elements of the splitting pattern
%  are below 2^53 (e.g. up to n = 55 for I = 1/2).

% see J.H.Freed, G.K.Fraenkel, J.Chem.Phys. 39(2), 326-348 (1963)
% https://doi.org/10.1063/1.1734250, eq. (4.32)

% Mathematical basis:
% Reduction of tensor product representations of rotation group
% using Clebsch-Gordan direct sum decomposition. For two spins:
% D^(j1)xD^(j2) = sum_{j=|j1-j2|}^{j1+j2} D^(j)
% For multiple spins, recursive.

function [F,N] = equivcouple(I,n)

if nargin==0, help(mfilename); return; end

if ~isnumeric(I) || ~isreal(I) || ~isscalar(I) || I<0 || mod(2*I,1)~=0
  error('I must be a nonnegative multiple of 1/2 (0, 1/2, 1, 3/2, ...).');
end
if ~isnumeric(n) || ~isreal(n) || ~isscalar(n) || n<0 || mod(n,1)~=0
  error('n must be a nonnegative integer (0, 1, 2, ...).');
end

% Special case n=0: a single spin-0
if n==0
  F = 0;
  N = 1;
  return
end

% Special case n=1: the single spin I is returned unchanged (the general
% code below would also return spins I-1, I-2, ... with zero abundance)
if n==1
  F = I;
  N = 1;
  return
end

% (1) List of reduced spin quantum numbers
F = I*n:-1:0;

% (2) Compute number of spins for each reduced spin quantum number:
% the number of spins F is the number of states with M=F minus the
% number of states with M=F+1 in the first-order splitting pattern
dPattern = diff([0, equivsplit(I,n)]);

N = dPattern(1:length(F));

end
