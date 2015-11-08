% This function computes the diffusion superoperator matrix in the LMKjK basis,
% for a diagonal diffusion tensor. It does not support ordering potentials.
% This function does not require a particular order of spatial basis functions.
%
%   Diff     principal values of diffusion tensor (s^-1)
%   basis    list of L,M,K,jK quantum numbers for spatial basis, one set per row
%   Gamma    diffusion superoperator in the given basis

function Gamma = diffsuperop(Diff,basis)

nBasis = size(basis,1);

L   = basis(:,1);
M   = basis(:,2);
K   = basis(:,3);
jK  = basis(:,4);

Rx = Diff(1);
Ry = Diff(2);
Rz = Diff(3);
Rd = 0.25*(Rx-Ry);
Rperp = 0.5*(Rx+Ry);

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
    
    if (K2==K1)
      val_ = Rperp*(L1*(L1+1)-K1^2) + Rz*K1^2;
    elseif (K2==K1+2)
      val_ = Rd*sqrt((1+K1zero)*(1+K2zero)*(L1+K2-1)*(L1+K2)*(L1-K2+1)*(L1-K2+2));
    elseif (K2==K1-2)
      val_ = Rd*sqrt((1+K1zero)*(1+K2zero)*(L1-K2-1)*(L1-K2)*(L1+K2+1)*(L1+K2+2));
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


