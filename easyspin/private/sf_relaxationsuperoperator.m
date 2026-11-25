% sf_relaxationsuperoperator  Relaxation superoperator for saffron (full engine)
%
%   [Gamma,minOverlap] = sf_relaxationsuperoperator(Sys,H)
%   [Gamma,minOverlap,GammaEigen,V] = sf_relaxationsuperoperator(Sys,H)
%
%   Builds the relaxation superoperator Gamma for the Liouville-space equation
%     d(sigma)/dt = -1i*[H,sigma] - Gamma*(sigma-sigma_eq)
%   where sigma is vectorized column-wise (sigma(:)).
%
%   Relaxation acts on the eigenstates of H. Relaxation times given as
%   matrices refer to product states, in the order of sop (descending m
%   for each spin). Each eigenstate is assigned to the product state it
%   overlaps most with (one-to-one).
%
%   Input:
%     Sys.Spins, Sys.nElectrons  spin quantum numbers (electrons first) and
%             number of electrons
%     Sys.T1  longitudinal relaxation time (us), scalar or nStates x nStates
%             matrix; T1(i,j) gives the rate 1/(2*T1(i,j)) from state j to
%             state i. For an isolated pair of states with symmetric rates,
%             the population difference relaxes with T1(i,j). If T1(i,j)
%             with i>j is zero, T1(j,i) is used. A scalar is used for all pairs.
%     Sys.T2  transverse relaxation time (us), scalar or nStates x nStates
%             matrix; T2(i,j) is the decay time of the coherence between
%             states i and j. Only the upper triangle is used. A scalar is
%             used for all electron coherences (between states that differ
%             in the state of at least one electron spin), but not for
%             nuclear coherences.
%             Zero means no relaxation, both for T1 and T2.
%     H       Hamiltonian (MHz) in the product basis, nStates x nStates
%
%   Output:
%     Gamma       relaxation superoperator (1/us) in the product basis
%     minOverlap  smallest squared overlap between an eigenstate and its
%                 assigned product state (1 for unmixed states)
%     GammaEigen  relaxation superoperator in the eigenbasis V, for the
%                 vectorized state V'*sigma*V
%     V           eigenvectors of H, from sf_eigenbasis

function [Gamma,minOverlap,GammaEigen,V] = sf_relaxationsuperoperator(Sys,H)

N = size(H,1);

% Population transfer rates W(i,j) from state j to state i
T1 = Sys.T1;
if isscalar(T1)
  T1 = T1*ones(N);
end
T1 = T1 - diag(diag(T1));
T1t = T1.';
useUpper = T1==0 & tril(true(N),-1);
T1(useUpper) = T1t(useUpper);
W = zeros(N);
W(T1>0) = 1./(2*T1(T1>0));

% Coherence decay rates
% A scalar T2 applies only to electron coherences, i.e. between product
% states that differ in the state of at least one electron spin. Nuclear
% coherences (same electron state) are not affected.
T2 = Sys.T2;
if isscalar(T2)
  nElectronStates = prod(2*Sys.Spins(1:Sys.nElectrons)+1);
  electronState = ceil((1:N)/(N/nElectronStates));  % electrons come first in Sys.Spins
  T2 = T2*(electronState.'~=electronState);
end
T2 = triu(T2,1);
T2 = T2 + T2.';
R2 = zeros(N);
R2(T2>0) = 1./T2(T2>0);

% Superoperator in the eigenbasis: rate equation for the populations
% (columns sum to zero, so the trace is conserved), exponential decay for
% the coherences
K = W - diag(sum(W,1));
Gamma = diag(R2(:));
idxPop = 1:N+1:N^2;
Gamma(idxPop,idxPop) = -K;

% Eigenstates of H, assigned one-to-one to product states
[V,~,minOverlap] = sf_eigenbasis(H);

% Transform to the product basis: vec(V'*sigma*V) = kron(V.',V')*vec(sigma)
GammaEigen = Gamma;
R = kron(V.',V');
Gamma = R'*GammaEigen*R;

end
