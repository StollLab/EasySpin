% This function computes the diffusion superoperator.
%
% This function does not require a particular order of spatial basis functions.

function Gamma = rdogamma(basis,DiffTensorValues,nSpin)

L   = basis(:,1);
M   = basis(:,2);
K   = basis(:,3);
jK  = basis(:,4);

Rx = DiffTensorValues(1);
Ry = DiffTensorValues(2);
Rz = DiffTensorValues(3);
Rd = 0.25*(Rx-Ry);
Rperp = 0.5*(Rx+Ry);

idx0 = 0;
for iBasis = 1:length(L)
  L1 = L(iBasis);
  M1 = M(iBasis);
  jK1 = jK(iBasis);
  K1 = K(iBasis);
  K1zero = (K1==0);
  idxr = (iBasis-1)*nSpin;
  
  for jBasis = iBasis:length(L)
    L2 = L(jBasis);
    if (L1~=L2), break; end
    M2 = M(jBasis);
    if (M1~=M2), continue; end
    jK2 = jK(jBasis);
    if (jK1~=jK2), continue; end
    K2 = K(jBasis);
    if (abs(K1-K2) > 2), continue; end
    
    K2zero = (K2==0);
    idxc = (iBasis-1)*nSpin;
    
    if (K2==K1)
      val = Rperp*(L1*(L1+1)-K1^2) + Rz*K1^2;
    elseif (K2==K1+2)
      val = Rd*sqrt((1+K1zero)*(1+K2zero)*(L1+K2-1)*(L1+K2)*(L1-K2+1)*(L1-K2+2));
    elseif (K2==K1-2)
      val = Rd*sqrt((1+K1zero)*(1+K2zero)*(L1-K2-1)*(L1-K2)*(L1+K2+1)*(L1+K2+2));
    else
      continue
    end
    
    idx_ = idx0 + (1:nSpin);
    bra(idx_) = idxr + (1:nSpin);
    ket(idx_) = idxc + (1:nSpin);
    el(idx_)  = val;
    
    idx0 = idx0 + nSpin;
    
  end
end

Gamma = sparse(bra,ket,el,length(L)*nSpin,length(L)*nSpin);
Gamma = triu(Gamma);

% Fill in lower triangular part
Gamma = Gamma + triu(Gamma,1).';

return
