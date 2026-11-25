% sf_eigenbasis  Eigenstates of a Hamiltonian, assigned to product states
%
%   [V,E,minOverlap] = sf_eigenbasis(H)
%
%   Diagonalizes H and assigns each eigenstate one-to-one to the product
%   state it overlaps most with. Column k of V is the eigenstate assigned
%   to product state k (in the order of sop, descending m for each spin).
%   Within sets of degenerate eigenstates, eig returns an arbitrary basis,
%   which is rotated to be closest to the assigned product states. Used by
%   the full engine of saffron, so that the relaxation superoperator
%   and the nuclear RF frame use the same eigenbasis.
%
%   Input:
%     H           Hamiltonian (MHz) in the product basis, nStates x nStates
%
%   Output:
%     V           eigenvectors (columns), in the order of the assigned
%                 product states
%     E           eigenvalues (MHz), same order, column
%     minOverlap  smallest squared overlap between an eigenstate and its
%                 assigned product state (1 for unmixed states)

function [V,E,minOverlap] = sf_eigenbasis(H)

N = size(H,1);

% Eigenstates of H, assigned one-to-one to product states by maximum overlap
H = (H+H')/2;
[V,E] = eig(H,'vector');
Overlap = abs(V).^2;  % Overlap(i,k): product state i, eigenstate k
pairs = matchpairs(1-Overlap,N);
idx = zeros(1,N);
idx(pairs(:,1)) = pairs(:,2);
V = V(:,idx);
E = E(idx);

% Within sets of degenerate eigenstates, rotate the basis to be closest to
% the assigned product states (polar decomposition of the overlap matrix)
tol = 1e-9*max(1,max(abs(E)));
[Es,order] = sort(E);
setBoundaries = [0; find(diff(Es)>tol); N];
for iSet = 1:numel(setBoundaries)-1
  k = order(setBoundaries(iSet)+1:setBoundaries(iSet+1));
  if numel(k)>1
    [Ul,~,Ur] = svd(V(k,k)');
    V(:,k) = V(:,k)*(Ul*Ur');
  end
end
minOverlap = min(abs(diag(V)).^2);

end
