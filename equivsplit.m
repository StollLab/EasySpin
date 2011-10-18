% equivsplit   Equivalent nuclei: EPR splitting pattern 
%
%  Ampl = equivsplit(I,n)
%
%  Computes the line intensity pattern of an EPR
%  spectrum due to an S=1/2 and n nuclear spins
%  with quantum number I. First-order perturbations
%  are assumed.
%
%  Example:
%
%     Ampl = equivsplit(1/2,5)
%
%  5 spins-1/2 give rise to a first-order splitting
%  pattern Ampl = [1 5 10 10 5 1] according to
%
%              1              0 spins
%            1   1            1 spin-1/2
%          1   2   1          2 spins-1/2
%        1   3   3   1        3 spins-1/2
%      1   4   6   4   1      4 spins-1/2
%    1   5  10  10   5   1    5 spins-1/2

function Ampl = equivsplit(I,n)

if (nargin==0), help(mfilename); return; end

% Some error checking
%-------------------------------------------------
if (nargin<2) || (nargin>2), error('Wrong number of input arguments!'); end

if (n<1) || mod(n,1)
  error('n must be an integer greater than 0.');
end

if (numel(I)~=1) || ~isreal(I) || mod(I,0.5) || (I<=0)
  error('I must be a positive multiple of 1/2.');
end

% Construct intensity pattern
%------------------------------------------------
% RowVec is the n-th row from a generalization of Pascal's
% triangle. The computation is iterative, since there are
% no closed expressions available.
RowVec = 1;
for iRow = 1:n
  RowVec = [RowVec zeros(1,2*I)];
  for iCol = length(RowVec):-1:1
    RowVec(iCol) = sum(RowVec(max(1,iCol-2*I):iCol));
  end
end

Ampl = RowVec;

return
