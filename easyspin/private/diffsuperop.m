% diffsuperop   Calculation of diffusion superoperator
%
% This function computes the diffusion superoperator matrix:
% - It can use either LML or LMKjK basis, depending on what is given in basis.
% - It requires a diagonal diffusion tensor.
% - It does not require any particular order of spatial basis functions.
%
% Input:
%   basis     Structure with quantum numbers for orientational
%             basis. If .jK is given, the K-symmetrized (L,M,K,jK) basis is
%             used, otherwise the (L,M,K) basis is used. All field must have
%             the same number of elements.
%     .L      L quantum numbers
%     .M      M quantum numbers
%     .K      K quantum numbers
%     .jK     jK quantum numbers (optional) If given, use LMKjK basis,
%             if missing, use LMK basis.
%   R         array with 3 principal values of diffusion tensor (s^-1)
%   Potential (optional) structure with ordering potential coefficients in
%             Potential.xlk, as returned by chili_xlk
%
% Output:
%   Gamma    diffusion superoperator in the given LMK or LMKjK basis (s^-1)

function Gamma = diffsuperop(basis,R,Potential)

usePotential = nargin==3 && isfield(Potential,'lambda') && any(Potential.lambda(:));

useSymmetrizedBasis = isfield(basis,'jK') && ~isempty(basis.jK) && any(basis.jK);

L = basis.L;
M = basis.M;
K = basis.K;
nBasis = numel(L);

if useSymmetrizedBasis
  jK = basis.jK;
else
  jK = zeros(nBasis,1);
end

switch numel(R)
  case 1
    % isotropic
    Rz = R;
    Rperp = R;
    Rd = 0;
  case 3
    % principal values
    Rz = R(3);
    Rd = (R(1)-R(2))/4;
    Rperp = (R(1)+R(2))/2;
  otherwise
    error('Three principal values of diffusion tensor required.');
end

% Treat the cases of isotropic and axial diffusion tensors in the absence of
% an ordering potential. In these cases, the diffusion operator matrix is
% diagonal in both the LMK and LMKjK basis.
if ~usePotential && Rd==0
  diagonal = Rperp*(L.*(L+1)-K.^2) + Rz*K.^2;
  Gamma = spdiags(diagonal,0,nBasis,nBasis);
  return
end

% Potential-independent part of diffusion operator
%-------------------------------------------------------------------------------
% This works for both the LMK basis and the K-symmetrized basis.
Np = @(L,K) sqrt((L*(L+1)-K*(K+1))*(L*(L+1)-(K+1)*(K+2)));
Nm = @(L,K) sqrt((L*(L+1)-K*(K-1))*(L*(L+1)-(K-1)*(K-2)));
idx = 0;
for b1 = 1:nBasis
  L1 = L(b1);
  M1 = M(b1);
  K1 = K(b1);
  jK1 = jK(b1);
  ph = jK1*(-1)^(L1+K1);
  for b2 = b1:nBasis % run only over upper triangular part (matrix is symmetric)
    L2 = L(b2);
    if L1~=L2, continue; end
    M2 = M(b2);
    if M1~=M2, continue; end
    jK2 = jK(b2);
    if jK1~=jK2, continue; end
    K2 = K(b2);
    if K1~=K2 && K1~=K2+2 && K1~=K2-2, continue; end
    
    % Calculate matrix element
    if K1==K2
      val_ = Rperp*L2*(L2+1) + (Rz-Rperp)*K2^2;
    else
      val_ = 0;
    end
    
    if Rd~=0
      if useSymmetrizedBasis
        val2_ = Np(L2,K2)*(K1==K2+2) + Nm(L2,K2)*((K1==K2-2)+ph*(-K1==K2-2));
        val_ = val_ + Rd*val2_/sqrt((1+(K1==0))*(1+(K2==0)));
      else
        val_ = val_ + Rd*(Np(L2,K2)*(K1==K2+2)+Nm(L2,K2)*(K1==K2-2));
      end
    end
    
    % Store non-zero value
    if val_==0, continue; end
    idx  = idx + 1;
    row(idx) = b1;
    col(idx) = b2;
    values(idx) = val_;
    if b1~=b2 % store value in lower triangular part
      idx  = idx + 1;
      row(idx) = b2;
      col(idx) = b1;
      values(idx) = val_;
    end
    
  end
end
Gamma = sparse(row,col,values,nBasis,nBasis);


% Potential-dependent part of diffusion operator
%-------------------------------------------------------------------------------
if ~usePotential, return; end

idx = 0;
XLK = Potential.xlk;
xLmax = size(XLK,1)-1;
for b1 = 1:nBasis
  L1  = L(b1);
  M1  = M(b1);
  K1  = K(b1);
  jK1 = jK(b1);
  
  % run only over upper triangular part
  for b2 = b1:nBasis
    L2  = L(b2);
    M2  = M(b2);
    K2  = K(b2);
    jK2 = jK(b2);
    
    if M1~=M2, continue; end
    
    prefactorL = sqrt((2*L1+1)*(2*L2+1));
    
    if useSymmetrizedBasis
      
      % calculate jK-dependent prefactor
      if K1==0 && K2==0
        prefactorjK = sqrt(jK1)'*sqrt(jK2)/4;
      elseif K1~=0 && K2~=0
        prefactorjK = sqrt(jK1)'*sqrt(jK2)/2;
      else
        prefactorjK = sqrt(jK1)'*sqrt(jK2)/sqrt(8);
      end
      prefactor = prefactorL*prefactorjK;
      
      val_ = 0;
      for Lx = abs(L1-L2):min(xLmax,L1+L2)
        idx_xL = Lx+1;
        if ~any(XLK(idx_xL,:)), continue; end
        
        % calculate M-dependent 3j-symbol
        jjjxM = wigner3j(L1,Lx,L2,M1,0,-M1);
        if jjjxM==0, continue; end
        
        % calculate K-dependent factors, 1st term (K2-K1)
        if abs(K2-K1)>Lx
          term1 = 0;
        else
          xlk_ = XLK(idx_xL,(K1-K2)+Lx+1);
          jjjxK_ = wigner3j(L1,Lx,L2,K1,K2-K1,-K2);
          sign_ = (-1)^(K1-M1) + jK1*jK2*(-1)^(Lx+K1-M1);
          term1 = sign_*xlk_*jjjxK_;
        end
        
        % calculate K-dependent factors, 2nd term (K1+K2)
        if abs(K1+K2)>Lx
          term2 = 0;
        else
          xlk_ = XLK(idx_xL,(K1+K2)+Lx+1);
          jjjxK_ = wigner3j(L1,Lx,L2,K1,-K2-K1,K2);
          sign_ = jK1*(-1)^(Lx+L2-M1) + jK2*(-1)^(L2+K1+K2-M1);
          term2 = sign_*xlk_*jjjxK_;
        end
        
        % combine terms
        val_ = val_ + prefactor * jjjxM * (term1 + term2);
        
      end
      
    else
      
      val_ = 0;
      for Lx = abs(L1-L2):min(xLmax,L1+L2)
        xlk_ = XLK(Lx+1,(K1-K2)+Lx+1);
        if xlk_==0, continue; end
        
        % calculate M-dependent 3j-symbol
        jjjxM = wigner3j(L1,Lx,L2,M1,0,-M1);
        if jjjxM==0, continue; end
        
        % calculate K-dependent 3j-symbol
        if abs(K2-K1)>Lx, continue; end
        jjjxK = wigner3j(L1,Lx,L2,K1,K2-K1,-K2);
        if jjjxK==0, continue; end
        
        % combine terms
        val_ = val_ + xlk_ * jjjxM * jjjxK;
      end
      val_ = (-1)^(K1-M1) * prefactorL * val_;
      
    end
    
    % Store non-zero value
    if val_==0, continue; end
    idx = idx + 1;
    rowp(idx) = b1;
    colp(idx) = b2;
    valp(idx)  = val_;
    if b1~=b2
      idx = idx + 1;
      rowp(idx) = b2;
      colp(idx) = b1;
      valp(idx)  = val_;
    end
    
  end
end

Gamma = Gamma + sparse(rowp,colp,valp,nBasis,nBasis);

return
