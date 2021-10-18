% pqorder returns an index array that allows reordering of the Liouville
% spin basis from the standard descending mm order to ascending pq order,
% as used by the Cornell programs (Freed).
%
% Input:
%  Spins   vector of spin quantum numbers
%
% Output:
%  idx    index vector to convert from descending mm order to ascending pq order
%  pq     list of (p,q) quantum numbers
%  mm     list of (m1,m2) quantum numbers
%
% If M is an operator matrix a with basis functions in descending mm order, then
% M(idx,idx) is the same operator with basis functions in ascending pq order.
%
% pqorder is used by chili.

function [idx,pq,mm] = pqorder(Spins)

nSpins = numel(Spins);
nStates = 2*Spins+1;

% Construct m quantum numbers for Hilbert spin basis
n1 = prod(nStates);
n2 = 1;
for iSpin = 1:nSpins
  n1 = n1/nStates(iSpin);
  m_ = repmat(Spins(iSpin):-1:-Spins(iSpin),n1,n2);
  m(:,iSpin) = m_(:);
  n2 = n2*nStates(iSpin);
end

% Construct pq and mm quantum numbers for Liouville spin basis
idx = 1:2;
for iSpin = 1:nSpins
  m1 = repmat(m(:,iSpin),1,size(m,1));
  m2 = m1.';
  m1 = m1(:);
  m2 = m2(:);
  mm(:,idx) = [m1 m2];
  pq(:,idx) = [m1-m2 m1+m2];
  idx = idx + 2;
end

[dummy,idx] = sortrows(pq);

return
