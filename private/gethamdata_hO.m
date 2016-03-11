% Hamiltonian diagonalization (sparse and full)
% used by resfreqs_matrix and resfields

% B       static field vector
% Sys       Spin System
% sparse       1 for sparse matrix calculation
% idxT    transition index (single number per transition)
% nLevels number of lowest levels to return


% 

% V     eigenvectors of Hamiltonian
% E     eigenvalues of Hamiltonian (energies)
% dEdB  derivatives of energy with respect to B
% dE    energy differences for the transitions listed in idxT

function [V,E,dEdB,dE] = gethamdata_hO(B,Sys,sparse,idxT,nLevels)
% Compute eigenvalues and eigenvectors of Hamiltonian
if sparse
  H = sham(Sys,B,'sparse');
  [V,E] = eigs(H,nLevels);
  E = diag(E).';
  [E,idx_] = sort(E);
  V = V(:,idx_);
  if nargout >2
    %calculate energies for numerical derivative
    nB = B./abs(B);
    dB = 0.01;
    H = sham(Sys,B+dB*nB,'sparse');
    E2 = eigs(H,nLevels);
    E2 = sort(E2);
  end
else
  H = sham(Sys,B);
  [V,E] = eig(H);
  E = diag(E).';
  if (nLevels<numel(E))
    E = E(1:nLevels);
    V = V(:,1:nLevels);
  end
  if nargout >2
    %calculate energies for numerical derivative
    nB = B./norm(B);
    dB = 0.01;
    H = sham(Sys,B+dB*nB);
    E2 = eig(H);
    if (nLevels<numel(E2))
      E2 = E(1:nLevels);
    end
  end
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
  dEdB = (E2.'-E)/dB;
end

if (nargout>3)
  % Compute transition frequencies Ev-Eu
  M = E(:);
  M = M(:,ones(1,length(E)));
  dE = M.' - M;
  dE = dE(idxT);
end
return
