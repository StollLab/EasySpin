% Hamiltonian diagonalization (sparse and dense)
% used by resfreqs_matrix and resfields
%
% Inputs:
%   B       static field
%   n0      direction of magnetic field
%   Sys     spin system
%   sparse  1 for sparse matrix calculation
%   idxT    transition index (single number per transition)
%   nLevels number of lowest levels to return
% 
% Outputs:
%   V       eigenvectors of Hamiltonian
%   E       eigenvalues of Hamiltonian (energies)
%   dEdB    derivatives of energy with respect to B
%   dE      energy differences for the transitions listed in idxT

function [V,E,dEdB,dE] = gethamdata_hO(B,n0, Sys,sparse,idxT,nLevels)

dB = 0.01;  % field step size for numerical derivative, mT

% Compute eigenvalues and eigenvectors of Hamiltonian
if sparse
  H = ham(Sys,B*n0,'sparse');
  [V,E] = eigs(H,nLevels);
  E = diag(E).';
  if sum(abs(imag(E)))>1e-6, error('Imaginary energies obtained! Please report!'); end
  E = real(E);
  [E,idx_] = sort(E);
  V = V(:,idx_);
  if nargout >2
    % Calculate energies for numerical derivative
    dB = 0.01;
    H2 = ham(Sys,(B+dB)*n0,'sparse');
    E2 = eigs(H2,nLevels);
    if sum(abs(imag(E2)))>1e-6, error('Imaginary energies obtained! Please report!'); end
    E2 = sort(real(E2));
  end
else
  H = ham(Sys,B*n0);
  [V,E] = eig(H);
  E = diag(E).';
  if sum(abs(imag(E)))>1e-6, error('Imaginary energies obtained! Please report!'); end
  [E,idx_] = sort(real(E));
  V = V(:,idx_);
  if nLevels<numel(E)
    E = E(1:nLevels);
    V = V(:,1:nLevels);
  end
  if nargout >2
    % Calculate energies for numerical derivative   
    H2= ham(Sys,(B+dB)*n0);
    E2 = eig(H2);
    if sum(abs(imag(E2)))>1e-6, error('Imaginary energies obtained! Please report!'); end
    E2 = sort(real(E2));
    if (nLevels<numel(E2))
      E2 = E2(1:nLevels);
    end
  end
end

% Compute more correct eigenvectors for zero-field degeneracies
% uses linear approximation for field dependence, give only reasonable
% results for 0th and 1th order terms in field much larger than higher
% orders.
if B==0
  if nargout > 2
    G = H2-H; 
  else
    H2= ham(Sys,(B+dB)*n0);
    G = H2-H;
  end
  dE = abs(diff(E)).';
  tol = 1e3*eps*max(dE);
  blk = cumsum([1; dE>tol]);
  GG = V'*G*V;
  GG = (GG+GG')/2; % important: symmetrise
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
  dEdB = (E2.'-E)/dB;
end

% Compute transition energies Ev-Eu
if nargout>3
  dE = -(E.' - E);  % such that upper triangle is positive
  dE = dE(idxT);
end

end
