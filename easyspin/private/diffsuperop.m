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

function Gamma = diffsuperop(R,basis)

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

if numel(R)==1
  R = R*ones(1,3);
end

Rx = R(1);
Ry = R(2);
Rz = R(3);
Rd = (Rx-Ry)/4;
Rperp = (Rx+Ry)/2;

% Treat the cases of isotropic and axial diffusion tensors, where the
% diffusion operator matrix is diagonal in LMK and LMKjK.
if Rd==0
  diagonal = Rperp*L.*(L+1) + (Rz-Rperp)*K.^2;
  Gamma = spdiags(diagonal,0,nBasis,nBasis);
  return
end

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
    
    if (L1~=L2), break; end
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

% Fill in lower triangular part
Gamma = Gamma + triu(Gamma,1).';

return
