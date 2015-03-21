function basisList = generatebasis(Basis,nSpin)

%%-------------------------------------------------------------------------
% Generate array of basis functions and their indices in full Liouvillian
%   
%   L  increasing from 0 to Lmax
%   M  decresing from L to -L
%   K  decreasing from L to 0
%   jK decreasing from +1 to -1 for nonzero K
%
%     L  M  K  jK
%%-------------------------------------------------------------------------

Leven = 0:2:Basis.evenLmax;
Lodd  = 1:2:Basis.oddLmax;
Llist = sort([Leven Lodd]);

Kmax = Basis.Kmax;
Mmax = Basis.Mmax;

iBasis = 1;
for L = Llist
  for M = min(L,Mmax):-1:max(-L,-Mmax)
    for K = min(L,Kmax):-1:0
      if (K~=0)
        basisList(iBasis,:)   = [L M K  1];
        basisList(iBasis+1,:) = [L M K -1];
        iBasis = iBasis + 2;
      else
        basisList(iBasis,:) = [L M K (-1)^L];
        iBasis = iBasis + 1;
      end
    end
  end
end

return
