% fullham takes as input the orientational basis size, the isotropic and
% anisotropic RBOs, and the spin basis, and gives as output the full
% spin/space Hamiltonian in Liouville space.
% The Wigner 3j-symbols have been precomupted and stored by jjjsymbol.


function [H0,H1,H2] = ksymham(basis,Q0,Q1,Q2,jjj0,jjj1,jjj2)


Lmax = max(basis(:,1));
L   = basis(:,1);
M   = basis(:,2);
K   = basis(:,3);
jK  = basis(:,4);
index = basis(:,5);
nSpace = sum((2*(0:Lmax)+1).^2);
nSpin = length(Q0);
nBasis = nSpace*nSpin;
nNonZero = length(basis);


%--------------------------------------------------------------------------
% compute the non-zero elements in the K-symmetrized Zeeman basis
%--------------------------------------------------------------------------


%j3idx = @(L,MK) L^2 + L - MK + 1;
%jjj0val = @(L1,MK1,L2,MK2) jjj0(j3idx(L1,MK1),j3idx(L2,MK2));
%jjj1val = @(L1,MK1,L2,MK2) jjj1(j3idx(L1,MK1),j3idx(L2,MK2));
%jjj2val = @(L1,MK1,L2,MK2) jjj2(j3idx(L1,MK1),j3idx(L2,MK2));


%--------------------------------------------------------------------------
% Rank 0
%--------------------------------------------------------------------------

i = 1;
for iBasis = 1:nNonZero
  L_  =  L(iBasis);
  M_  =  M(iBasis);
  K_  =  K(iBasis);
  idx = index(iBasis) - 1;
  
  NL = (2*L_+1);
  jjjM = jjj0(L_^2+L_-M_+1,L_^2+L_-M_+1);
  jjjK = jjj0(L_^2+L_-K_+1,L_^2+L_-K_+1);
  
  spinblock = (-1)^(K_-M_) * NL * jjjM * jjjK * Q0;
  
  [row,col,val] = find(spinblock);
  row = row + idx;
  col = col + idx;
  
  for j = 1:numel(row)
    braH0(i)  = row(j);
    ketH0(i)  = col(j); %#ok<*AGROW>
    elH0(i)   = val(j);
    i = i + 1;
  end
end


%--------------------------------------------------------------------------
% Rank 1
%--------------------------------------------------------------------------


i = 1;
for iBasis = 1:nNonZero
  L1_  =  L(iBasis);
  M1_  =  M(iBasis);
  K1_  =  K(iBasis);
  jK1_ = jK(iBasis);
  idx1 = index(iBasis) - 1;
  
  if K1_ == 0
    deltaK1 = 1;
  else
    deltaK1 = 0;
  end
  
  for jBasis = 1:nNonZero
    L2_  =  L(jBasis);
    if (abs(L1_-L2_) > 1), continue; end;
    M2_  =  M(jBasis);
    if (abs(M1_-M2_) > 1), continue; end;
    K2_  =  K(jBasis);
    if (abs(K1_-K2_) > 1), continue; end;
    jK2_ = jK(jBasis);
    
    idx2 = index(jBasis) - 1;
    NL = sqrt((2*L1_+1)*(2*L2_+1));
    jjjM = jjj1(L1_^2+L1_-M1_+1,L2_^2+L2_-M2_+1);
    jjjKa = jjj1(L1_^2+L1_-K1_+1,L2_^2+L2_-K2_+1);
    
    if K2_ == 0
      deltaK2 = 1;
    else
      deltaK2 = 0;
    end
    
    idxM  = 2-( M1_-M2_);
    idxKa = 2-( K1_-K2_);
    idxKb = 2-(-K1_+K2_);
    
    % first term
    
    prefactor = (1/(2*sqrt((1+deltaK1)*(1+deltaK2))))...
      * sqrt(jK1_)'*sqrt(jK2_) * NL * jjjM;
    sign1 = (-1)^(K1_-M1_);
    sign2 = (-1)^(K2_-M1_);
    
    spinblock = prefactor * jjjKa * (sign1*Q1{idxM,idxKa}...
      + jK1_*jK2_*sign2*Q1{idxM,idxKb});
    
    % second term
    
    if K1_+K2_ <= 1
      jjjKb = jjj1(L1_^2+L1_-K1_+1,L2_^2+L2_+K2_+1);
      idxKc = 2-(-K1_-K2_);
      idxKd = 2-(K1_+K2_);
      sign3 = (-1)^(L2_-M1_);
      sign4 = (-1)^(L2_+K1_+K2_-M1_);
      
      spinblock = spinblock + ...
        prefactor * jjjKb * (jK1_*sign3*Q1{idxM,idxKc}...
        + jK2_*sign4*Q1{idxM,idxKd});
    end
      
    [row,col,val] = find(spinblock);
    row = row + idx1;
    col = col + idx2;
    indices = i:i+numel(row)-1;
    
    braH1(indices) = row;
    ketH1(indices) = col;
    elH1(indices)  = val;
    
    i = i + numel(row);
    
  end
end


%--------------------------------------------------------------------------
% Rank 2
%--------------------------------------------------------------------------


i = 1;
for iBasis = 1:nNonZero
  L1_  =  L(iBasis);
  M1_  =  M(iBasis);
  K1_  =  K(iBasis);
  jK1_ = jK(iBasis);
  idx1 = index(iBasis) - 1;
  
  if K1_ == 0
    deltaK1 = 1;
  else
    deltaK1 = 0;
  end
  
  for jBasis = 1:nNonZero
    L2_  =  L(jBasis);
    if (abs(L1_-L2_) > 2), continue; end;
    M2_  =  M(jBasis);
    if (abs(M1_-M2_) > 2), continue; end;
    K2_  =  K(jBasis);
    if (abs(K1_-K2_) > 2), continue; end;
    jK2_ = jK(jBasis);
    
    idx2 = index(jBasis) - 1;
    NL = sqrt((2*L1_+1)*(2*L2_+1));
    jjjM = jjj2(L1_^2+L1_-M1_+1,L2_^2+L2_-M2_+1);
    jjjKa = jjj2(L1_^2+L1_-K1_+1,L2_^2+L2_-K2_+1);
    
    if K2_ == 0
      deltaK2 = 1;
    else
      deltaK2 = 0;
    end
    
    idxM  = 3-( M1_-M2_);
    idxKa = 3-( K1_-K2_);
    idxKb = 3-(-K1_+K2_);
    
    % first term
    
    prefactor = (1/(2*sqrt((1+deltaK1)*(1+deltaK2))))...
      * sqrt(jK1_)'*sqrt(jK2_) * NL * jjjM;
    sign1 = (-1)^(K1_-M1_);
    sign2 = (-1)^(K2_-M1_);
    
    spinblock = prefactor * jjjKa * (sign1*Q2{idxM,idxKa}...
      + jK1_*jK2_*sign2*Q2{idxM,idxKb});
    
    % second term
      
    if K1_+K2_ <= 2
      jjjKb = jjj2(L1_^2+L1_-K1_+1,L2_^2+L2_+K2_+1);
      idxKc = 3-(-K1_-K2_);
      idxKd = 3-(K1_+K2_);
      sign3 = (-1)^(L2_-M1_);
      sign4 = (-1)^(L2_+K1_+K2_-M1_);
      
      spinblock = spinblock + ...
        prefactor * jjjKb * (jK1_*sign3*Q2{idxM,idxKc}...
        + jK2_*sign4*Q2{idxM,idxKd});
    end
      
    [row,col,val] = find(spinblock);
    row = row + idx1;
    col = col + idx2;
    indices = i:i+numel(row)-1;
    
    braH2(indices) = row;
    ketH2(indices) = col;
    elH2(indices)  = val;
    
    i = i + numel(row);
    
  end
end


H0 = sparse(braH0,ketH0,elH0,nBasis,nBasis);
H1 = sparse(braH1,ketH1,elH1,nBasis,nBasis);
H1 = (H1 - diag(diag(H1))).';
H2 = sparse(braH2,ketH2,elH2,nBasis,nBasis);
H2 = (H2 - diag(diag(H2))).';

return


