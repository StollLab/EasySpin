% sf_spintilt  Tilt of the electron spin by a weak transverse field
%
%   dS = sf_spintilt(H,V,Sops,Szf)
%   [dS,dSsec] = sf_spintilt(H,V,Sops,Szf)
%
%   A weak field b (mT) adds b*V to the Hamiltonian H. To first order in b,
%   it changes the spin operators S_k projected onto each electron mS
%   manifold m, P_m*S_k*P_m, by b*dS_k. With the eigenstates |a> of H
%   (energies E_a), for a and a' in the same manifold,
%
%     <a|dS_k|a'> = sum over c in other manifolds of
%                   <a|S_k|c><c|V|a'>/(E_a'-E_c) + <a|V|c><c|S_k|a'>/(E_a-E_c)
%
%   This is the first-order perturbation of the eigenstates by b*V,
%   restricted to the admixture of other manifolds, i.e. first order in
%   1/nu_e. For a single electron level per manifold (S = 1/2), <a|dS|a> is
%   the derivative of the spin expectation value of level a with respect to
%   b. For S = 1/2 with the Zeeman Hamiltonian bmagn*B*zL.'*g*S, it is
%   mS*(g.'*n - q*(q.'*g.'*n))/(|g.'*zL|*B), with q = g.'*zL/|g.'*zL| and the
%   field direction n of b (the tilt of the quantization axis).
%
%   The RF field of an ENDOR pulse varies slowly compared with the electron
%   Larmor frequency, so the electron spin follows it adiabatically within
%   its manifold. Through the hyperfine coupling, the tilt dS gives the RF
%   field an additional effect on the nuclei (hyperfine enhancement).
%   Couplings between levels within a manifold (nuclear frequencies) are
%   not included, since they are not fast compared with the RF field.
%
%   Each eigenstate is assigned to the eigenspace of Szf (mS manifold) in
%   which it has the largest weight, as in sf_secularham.
%
%   The secular Hamiltonian of sf_secularham replaces the eigenstates of
%   each manifold by their projections onto the eigenspace of Szf,
%   orthonormalized. dSsec is dS in these states: it has the same matrix
%   elements between them as dS between the eigenstates, and no elements
%   between different eigenspaces of Szf. For S = 1/2, dS itself has
%   elements cos(theta)*<a|dS|a> on the diagonal of the Szf eigenbasis
%   (theta: angle between the quantization axis and Szf), and dSsec has
%   <a|dS|a> there.
%
%   Input:
%     H     Hamiltonian (MHz), field along the quantization axis of Szf
%     V     perturbation per unit field (MHz/mT), -n*mu for a field along n
%           with the magnetic moment operator mu
%     Sops  cell array of spin operators
%     Szf   total electron spin operator along the static field
%           All matrices are in the same basis.
%
%   Output:
%     dS     cell array of the changes of the operators in Sops (per mT), in
%            the basis of the input, Hermitian
%     dSsec  the same in the states of the secular Hamiltonian (see above),
%            in the basis of the input, Hermitian

function [dS,dSsec] = sf_spintilt(H,V,Sops,Szf)

[U,E] = eig(full(H+H')/2,'vector');

% mS manifold of each eigenstate
[Uz,mz] = eig(full(Szf+Szf')/2,'vector');
[~,~,mzIdx] = unique(round(2*mz)/2);
weights = (mzIdx==1:max(mzIdx)).'*abs(Uz'*U).^2;  % manifolds x eigenstates
[~,manifold] = max(weights,[],1);
sameManifold = manifold(:)==manifold(:).';

% D(a,c) = 1/(E_a-E_c) for states in different manifolds
dE = E - E.';
D = zeros(numel(E));
D(~sameManifold) = 1./dE(~sameManifold);

% Secular states (columns of T, one per eigenstate): projection of the
% eigenstates of each manifold onto the eigenspace of Szf, orthonormalized
% by the closest unitary (polar factor), as in sf_secularham
secular = nargout>1;
if secular
  T = zeros(numel(E));
  for m = 1:max(mzIdx)
    rows = mzIdx==m;
    cols = manifold==m;
    if nnz(rows)~=nnz(cols)
      error('The electron spin states cannot be assigned one-to-one to the mS manifolds.');
    end
    [L,~,R] = svd(Uz(:,rows)'*U(:,cols));
    T(:,cols) = Uz(:,rows)*(L*R');
  end
end

Ve = U'*V*U;
dS = cell(size(Sops));
dSsec = cell(size(Sops));
for k = 1:numel(Sops)
  Sk = U'*Sops{k}*U;
  dSeig = (Sk*(Ve.*D.') + (Ve.*D)*Sk).*sameManifold;  % eigenbasis
  dSk = U*dSeig*U';
  dS{k} = (dSk+dSk')/2;
  if secular
    dSk = T*dSeig*T';
    dSsec{k} = (dSk+dSk')/2;
  end
end

end
