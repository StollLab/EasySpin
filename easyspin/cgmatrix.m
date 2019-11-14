% cgmatrix     Transformation matrix between coupled and uncoupled representations
%
%   U2C = cgmatrix(S1,S2)
%
%   This function calculates the matrix representing the transformation from
%   uncoupled to coupled representation for a spin system consisting of two
%   spins with quantum numbers S1 and S2.
%
%   To transform vectors and matrices from the uncoupled to the coupled
%   representation, use
%      U2C = cgmatrix(S1,S2);
%      psi_c = U2C*psi_u;
%      H_c = U2C*H_u*U2C';
%
%   Input:
%      S1, S2: spin quantum numbers (1/2, 1, 3/2, 2, 5/2, etc)
%
%   Output:
%      U2C:    transformation matrix
%
%  Basis ordering convention - uncoupled basis:
%     top level: decreasing mS1 = S1,...,-S1
%     next level: decreasing mS2 = S2,...,-S2
%     |mS1,mS2> = 
%       |S1,S2>, |S1,S2-1>, ..., |S1-1,S2>, |S1-1,S2-1>, ..., |-S1,-S2>
%  Basis ordering convention - coupled basis
%     top level:  decreasing total spin St = S1+S2...|S1-S2|
%     next level: decreasing mSt = St...-St
%     |Stot,mStot> = 
%       |S1+S2,S1+S2>,|S1+S2,S1+S2-1>,...|S1+S2-1,S1+S2-1>,...|abs(S1-S2),-abs(S1-S2)>
%
%  cgmatrix() uses clebschgordan() to calculate the matrix elements of U2C.

function U2C = cgmatrix(S1,S2)

switch nargin
  case 0
    help(mfilename);
    return
  case 2
  otherwise
    error('Wrong number of input parameters. cgmatrix needs two inputs: S1 and S2.');
end

if numel(S1)~=1 || numel(S2)~=1
  error('S1 and S2 must be single numbers.');
end

if mod(S1*2,1) || mod(S2*2,1) || (S1<0.5) || (S2<0.5)
  error('S1 and S2 must be spin quantum numbers 1/2, 1, 3/2, etc.')
end

% preallocate transformation matrix
nStates = (2*S1+1)*(2*S2+1);
U2C = zeros(nStates);

% populate transformation matrix with Clebsch-Gordan coefficients
iu = 1; % index of uncoupled basis function
for mS1 = S1:-1:-S1
  for mS2 = S2:-1:-S2
    ic = 1; % index of coupled basis function
    for Stot = S1+S2:-1:abs(S1-S2)
      for mStot = Stot:-1:-Stot
        U2C(ic,iu) = clebschgordan(S1,S2,Stot,mS1,mS2,mStot);
        ic = ic + 1;
      end
    end
    iu = iu + 1;
  end
end

return
