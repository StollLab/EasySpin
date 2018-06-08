% Calculate matrix elements of Wigner functions
%
%    <L1,M1,K1|D^L_MK|L2,M2,K2>
%
% D = angmomops(L,M,K,Lx,Mx,Kx)
%
% Inputs:
%   L, K, M:  vectors of L1,M1,K1 and L2,M2,K2
%   Lx, Mx, Kx, vectors of L,M,K for the operator
%
% Outputs:
%   D: cell array of operator matrices, one for each
%      element in Lx,Mx,Kx

function DLMK = wignerops_LMK(L,M,K,Lx,Mx,Kx)

nOps = numel(Lx);
nBasis = numel(L);
for iOp = 1:nOps
  Lx_ = Lx(iOp);
  if Lx_==0
    DLMK{iOp} = speye(nBasis,nBasis);
    continue
  end
  
  D_ = zeros(nBasis,nBasis);
  
  for b1 = 1:nBasis
    L1 = L(b1);
    M1 = M(b1);
    K1 = K(b1);
    for b2 = 1:nBasis
      L2 = L(b2);
      M2 = M(b2);
      K2 = K(b2);
      dM = M1-M2;
      dK = K1-K2;
      
      % Apply 3j selection rules
      if Mx(iOp)~=dM, continue; end
      if Kx(iOp)~=dK, continue; end
      if abs(L1-Lx_)>L2 || L2>L1+Lx_, continue; end
      
      v = wigner3j([L1 Lx_ L2],[-M1 dM M2]);
      if v==0, continue; end
      v = v*wigner3j([L1 Lx_ L2],[-K1 dK K2]);
      if v==0, continue; end
      D_(b1,b2) = sqrt((2*L1+1)*(2*L2+1))*(-1)^(K1-M1)*v;
      
    end
  end
  
  DLMK{iOp} = sparse(D_);
  
end

if nOps==1
  DLMK = DLMK{1};
end
