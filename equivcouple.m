% equivcouple   Coupling of equivalent spins 
%
%  [K,N] = equivcouple(I,n)
%
%  The states due to n spins-I can be coupled to give
%  a set of independent spins. Their quantum numbers
%  are returned in K, their respective abundance
%  in N.
%
%  Example:
%
%     [K,N] = equivcouple(1/2,5)
%
%  5 spins-1/2 give rise to a first-order splitting
%  pattern [1 5 10 10 5 1] (see the function
%  equivsplit). This can be decomposed into one
%  spin-5/2, four spin-3/2 and five spin-1/2 according to
%
%          5  5          5 spins-1/2
%       4  4  4  4       4 spins-3/2
%    1  1  1  1  1  1    1 spin-5/2
%   ------------------
%    1  5 10 10  5  1    sum
%
%  so K = [2.5 1.5 0.5] and N = [1 4 5].

function [K,N] = equivcouple(I,n)

if (nargin==0), help(mfilename); return; end

RowVec = equivsplit(I,n);

% List of recoupled spin quantum numbers
%------------------------------------------------
largestSpin = (length(RowVec)-1)/2;
K = largestSpin:-1:0;

% Number of spins for each recoupled spin quantum number
% -------------------------------------------------------
nKSpins = length(K);
N = [1 diff(RowVec(1:nKSpins))];

if N(end)==0
  K(end) = [];
  N(end) = [];
end

return

% test code
[K,N] = equivcouple(1/2,5);
TestOK(1) = all(K==[5/2 3/2 1/2]) & all(N==[1 4 5]);

[K,N] = equivcouple(1,2);
TestOK(2) = all(K==[2 1 0]) & all(N==[1 1 1]);

TestOK
