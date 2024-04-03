% Hamiltonian diagonalization (sparse and dense) used by resfreqs_matrix and resfields
%
% Inputs:
%   B       static field magnitude
%   H0      field-independent part of Hamiltonian
%   muzL    magnetic dipole operator along lab z axis
%   idxT    transition index (single number per transition)
%   nLevels number of lowest levels to return
%
% Outputs:
%   V     eigenvectors of Hamiltonian
%   E     eigenvalues of Hamiltonian (energies)
%   dEdB  derivatives of energy with respect to B
%   dE    energy differences for the transitions listed in idxT
%
% The Hamiltonian is H = H0 - B*muzL

function [V,E,dEdB,dE] = gethamdata(B,H0,muzL,idxT,nLevels)

% Compute eigenvalues and eigenvectors of Hamiltonian
H = H0 - B*muzL;
if issparse(H)
  [V,E] = eigs(H,nLevels);
else
  [V,E] = eig(H);
end
E = diag(E).';
if ~issorted(E)
  [E,idx_] = sort(E);
  V = V(:,idx_);
end
if nLevels<numel(E)
  E = E(1:nLevels);
  V = V(:,1:nLevels);
end

% Compute correct eigenvectors for zero-field degeneracies
if B==0
  dE = abs(diff(E)).';
  tol = 1e3*eps*max(dE);
  blk = cumsum([1; dE>tol]);
  GG = V'*(-muzL)*V;
  GG = (GG+GG')/2; % important: symmetrize
  VV_ = [];
  for k = 1:max(blk)
    ix = find(blk==k);
    [v,~] = eig(GG(ix,ix));
    VV_ = blkdiag(VV_,v);
  end
  V = V*VV_;
end

% Compute first derivative of energies with respect to field, dE/dB
if nargout>2
  dEdB = real(diag(V'*(-muzL)*V)).';
end

% Compute transition energies Ev-Eu
if nargout>3
  dE = -(E.' - E);  % such that upper triangle is positive
  dE = dE(idxT);
end

end
