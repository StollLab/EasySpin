% cgmatrix     Transformation matrix between coupled and uncoupled representations
%
%   U2C = cgmatrix(S1,S2)
%   [U2C,Smtot] = cgmatrix(S1,S2)
%   [U2C,Smtot,m12] = cgmatrix(S1,S2)
%   ... = cgmatrix(S1,S2,Stot)
%
%   This function calculates the matrix representing the transformation from
%   uncoupled to coupled representation for a spin system consisting of two
%   spins with quantum numbers S1 and S2.
%
%   To transform vectors and matrices from the uncoupled to the coupled
%   representation, use
%      U2C = cgmatrix(S1,S2);
%      psi_c = U2C*psi_u;   % convert vector from uncoupled to coupled basis
%      H_c = U2C*H_u*U2C';  % convert matrix from uncoupled to coupled basis
%
%   The rows of U2C represent the coupled eigenvectors in the uncoupled
%   basis, and the colums of U2C represent the uncoupled eigenvectors in
%   the coupled basis.
%
%   Input:
%      S1, S2: spin quantum numbers (1/2, 1, 3/2, 2, 5/2, etc.)
%      Stot:   (optional) return only coupled states with total spin Stot
%
%   Output:
%      U2C:    transformation matrix with elements <Stot,mStot|mS1,mS2>,
%              where Stot is the total coupled spin and mStot its
%              associated magnetic quantum number
%      Smtot:  array of coupled-basis quantum numbers Stot and mStot for
%              each state of the coupled basis
%      m12:    array of coupled-basis quantum numbers mS1 and mS2 for
%              each state of the uncoupled basis
%
%   Basis ordering of the uncoupled basis:
%     top level: decreasing mS1 = S1,...,-S1
%     next level: decreasing mS2 = S2,...,-S2
%     |mS1,mS2> = 
%        |S1,S2>, |S1,S2-1>, ..., |S1-1,S2>, |S1-1,S2-1>, ..., |-S1,-S2>
%   Basis ordering of the coupled basis:
%     top level:  decreasing total spin Stot = S1+S2...|S1-S2|
%     next level: decreasing mStot = Stot...-Stot
%     |Stot,mStot> = 
%        |S1+S2,S1+S2>,|S1+S2,S1+S2-1>,...|S1+S2-1,S1+S2-1>,...|abs(S1-S2),-abs(S1-S2)>
%
%  cgmatrix() uses clebschgordan() to calculate the matrix elements of U2C.

function [U2C,Smtot,m12] = cgmatrix(S1,S2,Stot)

switch nargin
  case 0
    help(mfilename);
    return
  case 2
    Stot = [];
  case 3
  otherwise
    error('Wrong number of input parameters. cgmatrix needs two inputs: S1 and S2.');
end

if numel(S1)~=1 || numel(S2)~=1
  error('S1 and S2 must be single numbers.');
end

if mod(S1*2,1) || mod(S2*2,1) || S1<0.5 || S2<0.5
  error('S1 and S2 must be spin quantum numbers 1/2, 1, 3/2, etc.')
end

if ~isempty(Stot)
  if numel(Stot)~=1 || mod(Stot*2,1) || Stot>S1+S2 || Stot<abs(S1-S2)
    error('Invalid third input (Stot).')
  end
else
  Stot = S1+S2:-1:abs(S1-S2);
end

nStates = sum(2*Stot+1);

% Preallocate transformation matrix
U2C = zeros(nStates);
Smtot = zeros(nStates,2);
m12 = zeros(nStates,2);

% Populate transformation matrix with Clebsch-Gordan coefficients
iu = 1; % index of uncoupled basis function
for mS1 = S1:-1:-S1
  for mS2 = S2:-1:-S2
    m12(iu,:) = [mS1 mS2];
    ic = 1; % index of coupled basis function
    for Stot_ = Stot
      for mStot = Stot_:-1:-Stot_
        U2C(ic,iu) = clebschgordan(S1,S2,Stot_,mS1,mS2,mStot);
        if iu==1
          Smtot(ic,:) = [Stot_ mStot];
        end
        ic = ic + 1;
      end
    end
    iu = iu + 1;
  end
end

end
