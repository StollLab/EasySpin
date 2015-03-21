% Compute diffusion superoperator

function Gamma = rdogamma(basis,DiffTensorValues,nSpin)

L   = basis(:,1);
M   = basis(:,2);
K   = basis(:,3);
jK  = basis(:,4);
index = basis(:,5);
Lmax = max(L);
nSpatialBasis = length(L);

Rx = DiffTensorValues(1);
Ry = DiffTensorValues(2);
Rz = DiffTensorValues(3);
Rd = 0.25*(Rx-Ry);
Rperp = 0.5*(Rx+Ry);

idx0 = 0;
for iBasis = 1:nSpatialBasis
  L1 = L(iBasis);
  M1 = M(iBasis);
  jK1 = jK(iBasis);
  K1 = K(iBasis);
  deltaK1 = (K1==0);
  idxr = index(iBasis) - 1;
  
  for jBasis = iBasis:nSpatialBasis
    L2 = L(jBasis);
    if (L1~=L2), break; end
    L_ = L1;
    M2 = M(jBasis);
    if (M1~=M2), continue; end
    jK2 = jK(jBasis);
    if (jK1~=jK2), continue; end
    K2 = K(jBasis);
    if (abs(K1-K2) > 2), continue; end
    
    deltaK2 = (K2==0);
    idxc = index(jBasis) - 1;
    
    if (K2==K1)
      val = Rperp*(L1*(L1+1)-K1^2) + Rz*K1^2;
    elseif (K2==K1+2)
      val = Rd*sqrt((1+deltaK1)*(1+deltaK2)*(L_+K2-1)*(L_+K2)*(L_-K2+1)*(L_-K2+2));
    elseif (K2==K1-2)
      val = Rd*sqrt((1+deltaK1)*(1+deltaK2)*(L_-K2-1)*(L_-K2)*(L_+K2+1)*(L_+K2+2));
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

nSpace = sum((2*(0:Lmax)+1).^2);
nBasis = nSpace*nSpin;
Gamma = sparse(bra,ket,el,nBasis,nBasis);
Gamma = Gamma + triu(Gamma,1).';

return
