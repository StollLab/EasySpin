% This function computes the diffusion superoperator matrix in the LMK or LMKjK basis,
% for a diagonal diffusion tensor. It does not require any particular order of
% spatial basis functions. It does not support ordering potentials.
%
% Input:
%   R        array with 3 principal values of diffusion tensor (s^-1)
%   basis    list of (L,M,K) or (L,M,K,jK) quantum numbers for orientational
%            basis, one set per row
%
% Output:
%   Gamma    diffusion superoperator in the given LMK or LMKjK basis (s^-1)

function Gamma = diffsuperop(R,basis,XLK)

usePotential = nargin==3 && any(XLK(:)~=0);

nBasis = size(basis,1);

L = basis(:,1);
M = basis(:,2);
K = basis(:,3);

jKbasis = size(basis,2)==4;
if jKbasis
  jK = basis(:,4);
else
  jK = zeros(nBasis,1);
end

if usePotential && ~jKbasis
  error('Cannot use LMK basis with ordering potential.');
end

switch numel(R)
  case 1
    R = R*ones(1,3);
  case 3
    % principal values
  otherwise
end

Rx = R(1);
Ry = R(2);
Rz = R(3);
Rd = (Rx-Ry)/4;
Rperp = (Rx+Ry)/2;

% Treat the cases of isotropic and axial diffusion tensors, where the
% diffusion operator matrix is diagonal in LMK and LMKjK.
if ~usePotential && Rd==0
  diagonal = Rperp*L.*(L+1) + (Rz-Rperp)*K.^2;
  Gamma = spdiags(diagonal,0,nBasis,nBasis);
  return
end

% Potential-independent part of diffusion operator
%-------------------------------------------------------------------------------
idx = 1;
for b1 = 1:nBasis
  L1  = L(b1);
  M1  = M(b1);
  K1  = K(b1);
  jK1 = jK(b1);
  
  K1zero = (K1==0);
  
  for b2 = b1:nBasis
    L2  = L(b2);
    M2  = M(b2);
    K2  = K(b2);
    jK2 = jK(b2);
    
    if (L1~=L2), continue; end
    if (M1~=M2), continue; end
    if (abs(K1-K2) > 2), continue; end
    
    if (jK1~=jK2), continue; end
    
    K2zero = (K2==0);
    
    if (K1==K2)
      val_ = Rperp*(L1*(L1+1)-K1^2) + Rz*K1^2;
    elseif (K1==K2+2)
      if jKbasis
        val_ = Rd*sqrt((1+K1zero)*(1+K2zero)*(L1-K2-1)*(L1-K2)*(L1+K2+1)*(L1+K2+2));
      else
        val_ = Rd*sqrt((L1-K2-1)*(L1-K2)*(L1+K2+1)*(L1+K2+2));
      end
    elseif (K1==K2-2)
      if jKbasis
        val_ = Rd*sqrt((1+K1zero)*(1+K2zero)*(L1+K2-1)*(L1+K2)*(L1-K2+1)*(L1-K2+2));
      else
        val_ = Rd*sqrt((L1+K2-1)*(L1+K2)*(L1-K2+1)*(L1-K2+2));
      end
    else
      continue
    end
    
    bra(idx) = b1;
    ket(idx) = b2;
    val(idx)  = val_;
    idx = idx + 1;
    
  end
end
Gamma = sparse(bra,ket,val,nBasis,nBasis);

% Potential-dependent part of diffusion operator
%-------------------------------------------------------------------------------
if usePotential
  idx = 1;
  xLmax = size(XLK,1)-1;
  for b1 = 1:nBasis
    L1  = L(b1);
    M1  = M(b1);
    K1  = K(b1);
    jK1 = jK(b1);
    
    for b2 = b1:nBasis
      L2  = L(b2);
      M2  = M(b2);
      K2  = K(b2);
      jK2 = jK(b2);
      
      % calculate prefactors
      if (M1~=M2)
        prefactorL = 0;
      else
        prefactorL = sqrt((2*L1+1)*(2*L2+1));
      end
      
      if (K1==0 && K2==0)
        prefactorK = sqrt(jK1)'*sqrt(jK2)/4;
      elseif (K1~=0 && K2~=0)
        prefactorK = sqrt(jK1)'*sqrt(jK2)/2;
      else
        prefactorK = sqrt(jK1)'*sqrt(jK2)/sqrt(8);
      end
      
      val_ = 0;
      prefactor = prefactorL*prefactorK;
      for xL = abs(L1-L2):min(xLmax,L1+L2)
        idx_xL = xL+1;
        if ~any(XLK(idx_xL,:)), continue; end
        
        % calculate M-dependent 3j-symbol
        if abs(M1)>L1 || abs(M1)>L2
          jjjxM = 0;
        else
          jjjxM = wigner3j(L1,xL,L2,M1,0,-M1);
        end
        
        % 1st term (K2-K1)
        if abs(K1-K2)>xL
          xlk_1 = 0;
        else
          idx_xK_1 = (K1-K2)+xL+1;
          xlk_1 = XLK(idx_xL,idx_xK_1);
        end
        if abs(K1)>L1 || abs(K2-K1)>xL || abs(K2)>L2
          jjjxK_1 = 0;
        else
          jjjxK_1 = wigner3j(L1,xL,L2,K1,K2-K1,-K2);
        end
        sign1 = (-1)^(K1-M1) + jK1*jK2*(-1)^(xL+K1-M1);
        
        % 2nd term (K2-K1)
        if abs(K1+K2)>xL
          xlk_2 = 0;
        else
          idx_xK_2 = (K1+K2)+xL+1;
          xlk_2 = XLK(idx_xL,idx_xK_2);
        end
        if abs(K1)>L1 || abs(K2+K1)>xL || abs(K2)>L2
          jjjxK_2 = 0;
        else
          jjjxK_2 = wigner3j(L1,xL,L2,K1,-K2-K1,K2);
        end
        sign2 = jK1*(-1)^(xL+L2-M1) + jK2*(-1)^(L2+K1+K2-M1);
        
        % combine terms
        val_ = val_ + prefactor * jjjxM * ...
               (sign1*xlk_1*jjjxK_1 + sign2*xlk_2*jjjxK_2);
      end
      
      bra(idx) = b1;
      ket(idx) = b2;
      val(idx)  = val_;
      idx = idx + 1;
    end
  end
  Gamma_U = sparse(bra,ket,val,nBasis,nBasis);
  Gamma = Gamma + Gamma_U;
end

% Fill in lower triangular part
Gamma = Gamma + triu(Gamma,1).';

return
