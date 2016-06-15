% equivsplit   Equivalent nuclei: EPR splitting pattern 
%
%  Ampl = equivsplit(I,n)
%
%  Computes the line intensity pattern of an EPR
%  spectrum due to an S=1/2 and n equivalent nuclear spins
%  with quantum number I. First-order perturbations
%  are assumed.
%
%  Example:
%
%     Intensity = equivsplit(1/2,5)
%
%  5 spins-1/2 give rise to a first-order splitting
%  pattern Intensity = [1 5 10 10 5 1] according to
%
%              1              0 spins
%            1   1            1 spin-1/2
%          1   2   1          2 spins-1/2
%        1   3   3   1        3 spins-1/2
%      1   4   6   4   1      4 spins-1/2
%    1   5  10  10   5   1    5 spins-1/2

function Intensity = equivsplit(I,n)

if (nargin==0), help(mfilename); return; end

% Some error checking
%-------------------------------------------------
if (nargin<2) || (nargin>2), error('Wrong number of input arguments!'); end

if (n<1) || mod(n,1)
  error('The second input argument, n, must be an integer greater than 0.');
end

if (numel(I)~=1) || ~isreal(I) || mod(I,0.5) || (I<0)
  error('The spin quantum number I (first input argument) must be a positive multiple of 1/2.');
end

% Construct amplitude pattern
%-------------------------------------------------------------------------
% Intensity is the n-th row from a generalization of Pascal's triangle.
% The computation is brute-force iterative, as I couldn't find closed
% expressions for the generalized binomial coefficients of order 2*I+1 that
% constitute the (n+1)-th row of the associated the generalized Pascal triangle.

E = ones(1,2*I+1);
Intensity = 1;
for q = 1:n
  Intensity = conv(Intensity,E);
end

return
