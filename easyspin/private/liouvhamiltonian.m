% This function takes as input the orientational basis size, the isotropic and
% anisotropic RBOs, and the spin basis, and gives as output the full
% spin/space Hamiltonian in Liouville space.
% The Wigner 3j-symbols have been precomupted by jjjsymbol and stored in
% jjj0, jjj1, and jjj2.

% This function does not require a particular order of spatial basis functions.

function H = liouvhamiltonian(basis,Q0,Q1,Q2,jjj0,jjj1,jjj2)

L   = basis.L;
M   = basis.M;
K   = basis.K;
jK  = basis.jK;
nSpinBasis = length(Q0);
nOriBasis = numel(L);
nTotalBasis = nOriBasis*nSpinBasis;

%--------------------------------------------------------------------------
% Rank 0
%--------------------------------------------------------------------------
[row,col,val] = find(Q0);
indices = 1:numel(val);
for iBasis = 1:nOriBasis
  L_  =  L(iBasis);
  M_  =  M(iBasis);
  K_  =  K(iBasis);
  idx0 = (iBasis-1)*nSpinBasis;
  
  jjjM = jjj0(L_^2+L_+M_+1,L_^2+L_-M_+1);
  jjjK = jjj0(L_^2+L_+K_+1,L_^2+L_-K_+1);
  
  elH0(indices)  = (-1)^(K_-M_) * (2*L_+1) * jjjM * jjjK * val;
  braH0(indices) = idx0 + row;
  ketH0(indices) = idx0 + col;
  
  indices = indices + numel(val);
end

% Assemble sparse matrix
H0 = sparse(braH0,ketH0,elH0,nTotalBasis,nTotalBasis);


%--------------------------------------------------------------------------
% Rank 1
%--------------------------------------------------------------------------
if ~isempty(Q1)
  i = 1;
  for iBasis = 1:nOriBasis
    L1 = L(iBasis);
    M1 = M(iBasis);
    K1 = K(iBasis);
    jK1 = jK(iBasis);
    idx1 = (iBasis-1)*nSpinBasis;
    
    K1zero = (K1==0);
    
    for jBasis = 1:nOriBasis
      L2 = L(jBasis);
      if (abs(L1-L2) > 1), continue; end
      M2 = M(jBasis);
      if (abs(M1-M2) > 1), continue; end
      K2 = K(jBasis);
      if (abs(K1-K2) > 1), continue; end
      jK2 = jK(jBasis);
      
      idx2 = (jBasis-1)*nSpinBasis;
      NL = sqrt((2*L1+1)*(2*L2+1));
      jjjM = jjj1(L1^2+L1-M1+1,L2^2+L2-M2+1);
      jjjKa = jjj1(L1^2+L1-K1+1,L2^2+L2-K2+1);
      
      K2zero = (K2==0);
      
      idxM  = 2-( M1-M2);
      idxKa = 2-( K1-K2);
      idxKb = 2-(-K1+K2);
      
      % first term
      prefactor = (1/(2*sqrt((1+K1zero)*(1+K2zero))))...
        * sqrt(jK1)'*sqrt(jK2) * NL * jjjM;
      
      spinblock = prefactor * jjjKa * ...
        ((-1)^(K1-M1)*Q1{idxM,idxKa} + jK1*jK2*(-1)^(K2-M1)*Q1{idxM,idxKb});
      
      % second term
      if (K1+K2<=1)
        jjjKb = jjj1(L1^2+L1-K1+1,L2^2+L2+K2+1);
        idxKc = 2-(-K1-K2);
        idxKd = 2-(K1+K2);
        
        spinblock = spinblock + ...
          prefactor * jjjKb * (jK1*(-1)^(L2-M1)*Q1{idxM,idxKc}...
          + jK2*(-1)^(L2+K1+K2-M1)*Q1{idxM,idxKd});
      end
      
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
  
  % Fill in lower triangular part
  H1 = H1 + triu(H1,1).';
else
  H1 = sparse(0);
end

%--------------------------------------------------------------------------
% Rank 2
%--------------------------------------------------------------------------

i = 1;
for iBasis = 1:nOriBasis
  L1 = L(iBasis);
  M1 = M(iBasis);
  K1 = K(iBasis);
  jK1 = jK(iBasis);
  idx1 = (iBasis-1)*nSpinBasis;
  
  K1zero = (K1==0);
  
  for jBasis = iBasis:nOriBasis
    L2 = L(jBasis); % L2 >= L1 if basis set is L-ordered
    if (abs(L1-L2)>2), break; end
    M2 = M(jBasis);
    if (abs(M1-M2)>2), continue; end
    K2 = K(jBasis);
    if (abs(K1-K2)>2), continue; end
    jK2 = jK(jBasis);
    
    idx2 = (jBasis-1)*nSpinBasis;
    
    NL = sqrt((2*L1+1)*(2*L2+1));
    jjjM = jjj2(L1^2+L1-M1+1,L2^2+L2+M2+1);
    jjjKa = jjj2(L1^2+L1-K1+1,L2^2+L2+K2+1);
    
    K2zero = (K2==0);
    
    idxM  = 3-( M1-M2);
    idxKa = 3-( K1-K2);
    idxKb = 3-(-K1+K2);
    
    prefactor = (1/(2*sqrt((1+K1zero)*(1+K2zero))))...
      * sqrt(jK1)'*sqrt(jK2) * NL * jjjM;
    
    % first term
    spinblock = prefactor * jjjKa * ...
      ((-1)^(K1-M1)*Q2{idxM,idxKa} + jK1*jK2*(-1)^(K2-M1)*Q2{idxM,idxKb});
    
    % second term
    if (K1+K2<=2)
      jjjKb = jjj2(L1^2+L1-K1+1,L2^2+L2-K2+1);
      idxKc = 3-(-K1-K2);
      idxKd = 3-(K1+K2);
      
      spinblock = spinblock + ...
        prefactor * jjjKb * (jK1*(-1)^(L2-M1)*Q2{idxM,idxKc}...
        + jK2*(-1)^(L2+K1+K2-M1)*Q2{idxM,idxKd});
    end
    
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
H2 = H2 + triu(H2,1).';

% Combine rank-0, rank-1 and rank-2 parts of Hamiltonian
H = H0 + H1 + H2;

return
