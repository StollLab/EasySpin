

function basis = generatebasis(LLKM,nSpin)


%%-------------------------------------------------------------------------
% Generate array of basis functions and their indices in full Liouvillian
%   
%   L  increasing from 0 to Lmax
%   M  decresing from L to -L
%   K  decreasing from L to 0
%   jK decreasing from +1 to -1 for nonzero K
%
%     L  M  K  jK  index
%%-------------------------------------------------------------------------


if mod(LLKM(1),2) ~= 0
  LLKM(1) = LLKM(1) + 1;
end
if mod(LLKM(2),2) == 0
  LLKM(2) = LLKM(2) + 1;
end

Lmax  = max(LLKM(1),LLKM(2));
Leven = 0:2:LLKM(1);
Lodd  = 1:2:LLKM(2);
L = sort([Leven Lodd]);

Kmax = LLKM(3);
Mmax = LLKM(4);

iBasis = 1;
for L_ = L
  idxL = sum((2*(0:L_-1)+1).^2) * nSpin;
  for M_ = min(L_,Mmax):-1:max(-L_,-Mmax)
    idxM = (L_-M_)*(2*L_+1) * nSpin;
    for K_ = min(L_,Kmax):-1:0
      idxK = 2*(L_-K_) * nSpin;
      idx = idxL + idxM + idxK + 1;
      if K_ ~= 0
        basis(iBasis,:)   = [L_ M_ K_  1 idx];
        basis(iBasis+1,:) = [L_ M_ K_ -1 idx+nSpin];
        iBasis = iBasis + 2;
      else
        jK_ = (-1)^L_;
        basis(iBasis,:) = [L_ M_ K_ jK_ idx];
        iBasis = iBasis + 1;
      end
    end
  end
end


return