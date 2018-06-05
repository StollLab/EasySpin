% This function takes as input the orientational basis size, the isotropic and
% anisotropic RBOs, and the spin basis, and gives as output the full
% spin/space Hamiltonian in Liouville space.
% The Wigner 3j-symbols have been precomupted by jjjsymbol and stored in
% jjj0, jjj1, and jjj2.

% This function does not require a particular order of spatial basis functions.

function H = liouvhamiltonian_LMK(basis,Q0,Q1,Q2,jjj0,jjj1,jjj2)

if isfield(basis,'jK') && ~isempty(basis.jK)
  error('This function expects an LMK basis, without jK.');
end

L = basis.L;
M = basis.M;
K = basis.K;
nSpinBasis = length(Q0);
nOriBasis = numel(L);
nTotalBasis = nOriBasis*nSpinBasis;

%--------------------------------------------------------------------------
% Rank 0
%--------------------------------------------------------------------------
% only one loop, since the matrix is diagonal (L1=L2=L, M1=M2=M, K1=K2=K)
[row,col,Qval] = find(Q0);
indices = 1:numel(Qval);
for b1 = 1:nOriBasis
  L_  =  L(b1);
  M_  =  M(b1);
  K_  =  K(b1);
  idx0 = (b1-1)*nSpinBasis;
  
  jjjM = jjj0(L_^2+L_+M_+1,L_^2+L_-M_+1); % (L 0 L; M 0 -M)
  jjjK = jjj0(L_^2+L_+K_+1,L_^2+L_-K_+1); % (L 0 L; K 0 -K)
  
  elH0(indices)  = (-1)^(K_-M_) * (2*L_+1) * jjjM * jjjK * Qval;
  braH0(indices) = idx0 + row;
  ketH0(indices) = idx0 + col;
  
  indices = indices + numel(Qval);
end

% Assemble sparse matrix
H0 = sparse(braH0,ketH0,elH0,nTotalBasis,nTotalBasis);


%--------------------------------------------------------------------------
% Rank 1
%--------------------------------------------------------------------------
if ~isempty(Q1)
  i = 1;
  for b1 = 1:nOriBasis
    L1 = L(b1);
    M1 = M(b1);
    K1 = K(b1);
    idx1 = (b1-1)*nSpinBasis;
    
    for b2 = b1:nOriBasis
      L2 = L(b2);
      if abs(L1-L2) > 1, continue; end
      M2 = M(b2);
      if abs(M1-M2) > 1, continue; end
      K2 = K(b2);
      if abs(K1-K2) > 1, continue; end
      
      idx2 = (b2-1)*nSpinBasis;
      NL = sqrt((2*L1+1)*(2*L2+1));
      jjjM = jjj1(L1^2+L1-M1+1,L2^2+L2-M2+1); % (L1 1 L2; M1 M2-M1 -M1)
      jjjK = jjj1(L1^2+L1-K1+1,L2^2+L2-K2+1); % (L1 1 L2; K1 K2-K1 -K1)
      idxM = 2-(M1-M2);
      idxK = 2-(K1-K2);
      spinblock = (-1)^(K1-M1) * NL * jjjM * jjjK * Q1{idxM,idxK};
      
      % Store values
      [row,col,val] = find(spinblock);
      indices = i:i+numel(row)-1;
      braH1(indices) = row + idx1;
      ketH1(indices) = col + idx2;
      elH1(indices)  = val;
      i = i + numel(row);
      
    end
  end
  % Assemble sparse matrix
  H1 = sparse(braH1,ketH1,elH1,nTotalBasis,nTotalBasis);
  
  % Set all elements below the diagonal to zero
  H1 = triu(H1);
  
  % Fill in lower triangular part
  H1 = H1 + triu(H1,1).';
else
  H1 = sparse(0);
end

%--------------------------------------------------------------------------
% Rank 2
%--------------------------------------------------------------------------
i = 1;
for b1 = 1:nOriBasis
  L1 = L(b1);
  M1 = M(b1);
  K1 = K(b1);
  idx1 = (b1-1)*nSpinBasis;
  
  for b2 = b1:nOriBasis
    L2 = L(b2);
    if abs(L1-L2)>2, continue; end
    M2 = M(b2);
    if abs(M1-M2)>2, continue; end
    K2 = K(b2);
    if abs(K1-K2)>2, continue; end
    
    idx2 = (b2-1)*nSpinBasis;
    NL = sqrt((2*L1+1)*(2*L2+1));
    jjjM = jjj2(L1^2+L1-M1+1,L2^2+L2+M2+1); % (L1 2 L2; M1 M2-M1 -M2)
    jjjK = jjj2(L1^2+L1-K1+1,L2^2+L2+K2+1); % (L1 2 L2; K1 K2-K1 -K2)
    idxM = 3-(M1-M2);
    idxK = 3-(K1-K2);
    spinblock = (-1)^(K1-M1) * NL * jjjM * jjjK * Q2{idxM,idxK};
    
    % Store values
    [row,col,val] = find(spinblock);
    indices = i:i+numel(val)-1;
    braH2(indices) = row + idx1;
    ketH2(indices) = col + idx2;
    elH2(indices)  = val;
    i = i + numel(val);
    
  end
end

% Assemble sparse matrix
H2 = sparse(braH2,ketH2,elH2,nTotalBasis,nTotalBasis);

% Set all elements below the diagonal to zero
H2 = triu(H2);

% Fill in lower triangular part
H2 = H2 + triu(H2,1)';

% Combine rank-0, rank-1 and rank-2 parts of Hamiltonian
H = H0 + H1 + H2;

return
