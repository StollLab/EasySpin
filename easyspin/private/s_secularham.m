% s_secularham  Decouple electron mS manifolds of a lab-frame Hamiltonian
%
%   [Hsec,minWeight] = s_secularham(H,Spins,nElectrons)
%
%   Input:
%     H           Hamiltonian (MHz) in the lab frame, field along z
%     Spins       spin quantum numbers (Sys.Spins), electrons listed first
%     nElectrons  number of electron spins
%
%   Output:
%     Hsec        Hamiltonian (MHz) that commutes with total electron Sz
%     minWeight   for each mS manifold (in increasing order of mS), the
%                 smallest weight of an assigned eigenvector within the
%                 manifold (1 = no mixing with other manifolds)
%
% Returns a Hamiltonian with the same eigenvalues as H but no matrix
% elements between different total electron mS values, so that it can be
% transformed into a frame rotating about z. The eigenvectors of H are
% projected onto their dominant mS manifold and then minimally adjusted
% to be orthonormal again. Energy shifts caused by the removed couplings
% are retained.
%
% Requires the electron Zeeman interaction to dominate. Energies are kept
% exactly, but admixtures between manifolds (e.g. from g anisotropy, which
% tilts the quantization axis away from z) are removed from the
% eigenvectors, which changes transition intensities slightly.
%
% An error is issued if eigenvectors cannot be assigned one-to-one to the
% manifolds, or if any weight is below 0.7, since then the projection
% changes the eigenvectors too much. Callers should warn if a weight is
% below 0.9 (admixture amplitude above about 0.3).

function [Hsec,minWeight] = s_secularham(H,Spins,nElectrons)

% Construct total electron Sz and get total mS for each state
Sz_total = 0;
for iSpin = 1:nElectrons
  Sz_total = Sz_total + sop(Spins,[iSpin,3]);
end
mS_total = diag(Sz_total);

% Manifold index via total electron mS
[~,~,manifold] = unique(mS_total);
nManifolds = max(manifold);
inManifold = manifold==1:nManifolds;  % nStates x nManifolds indicator matrix

% Eigenvectors and energies of the full lab-frame Hamiltonian
[V,E] = eig((H+H')/2,'vector');

% Calculate weights of all eigenvectors in all mS manifolds. Assign each
% eigenvector to the manifold it has the largest weight in. Error if this
% assignment is not one-to-one, or if any assigned weight is too small.
[weight,assigned] = max(inManifold.'*abs(V).^2,[],1);
separable = isequal(sort(assigned(:)),sort(manifold));
minWeight = accumarray(assigned(:),weight(:),[nManifolds 1],@min);
if ~separable || any(minWeight<0.7)
  error(['The electron spin states are too strongly mixed to be separated ' ...
    'into mS manifolds (the electron Zeeman interaction does not dominate). ' ...
    'Use Opt.SimFreq = 0 to simulate in the lab frame.']);
end

% For each manifold, project assigned eigenvectors into manifold,
% re-orthonormalize them via SVD, and rebuild Hamiltonian block using
% exact energies.
Hsec = zeros(size(H));
for m = 1:nManifolds
  rows = manifold==m;
  cols = assigned==m;
  [L,~,R] = svd(V(rows,cols));
  U = L*R';
  Hsec(rows,rows) = U*diag(E(cols))*U';
end
Hsec = (Hsec+Hsec')/2;

end
