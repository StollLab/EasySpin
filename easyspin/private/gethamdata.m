% Hamiltonian diagonalization (sparse and full)
% used by resfreqs_matrix and resfields

% B       static field magnitude
% F       field-independent part of Hamiltonian
% G       field-dependent part of Hamiltonian
% idxT    transition index (single number per transition)
% nLevels number of lowest levels to return


% The Hamiltonian is H = F + B*G

% V     eigenvectors of Hamiltonian
% E     eigenvalues of Hamiltonian (energies)
% dEdB  derivatives of energy with respect to B
% dE    energy differences for the transitions listed in idxT

function [V,E,dEdB,dE] = gethamdata(B,F,G,idxT,nLevels)

% Compute eigenvalues and eigenvectors of Hamiltonian
if issparse(F)
  [V,E] = eigs(F+B*G,nLevels);
else
  [V,E] = eig(F+B*G);
end
E = diag(E).';
if ~issorted(E)
  [E,idx_] = sort(E);
  V = V(:,idx_);
end
if (nLevels<numel(E))
  E = E(1:nLevels);
  V = V(:,1:nLevels);
end
% Compute correct eigenvectors for zero-field degeneracies
if (B==0)
  dE = abs(diff(E)).';
  tol = 1e3*eps*max(dE);
  blk = cumsum([1; dE>tol]);
  GG = V'*G*V;
  GG = (GG+GG')/2; % important: symmetrise
  VV_ = [];
  for k=1:max(blk)
    ix = find(blk==k);
    [v,e] = eig(GG(ix,ix));
    VV_ = blkdiag(VV_,v);
  end
  V = V*VV_;
end

if (nargout>2)
  % Compute first derivative of energies with respect to field, dE/dB
  dEdB = real(diag(V'*G*V)).';
end

if (nargout>3)
  % Compute transition frequencies Ev-Eu
  M = repmat(E(:),1,numel(E));
  dE = M.' - M;
  dE = dE(idxT);
end

return
