function [nRows,Indices] = chili_basiscount(Sys,Basis,DirTilt)

evenLmax = Basis.evenLmax;
oddLmax = Basis.oddLmax;
deltaL = Basis.deltaL;

jkmn = Basis.jKmin;
Kmax = Basis.Kmax;
deltaK = Basis.deltaK;

Mmax = Basis.Mmax;

pSmin = Basis.pSmin;
pImax = Basis.pImax;

I = Sys.I;

iRow = 0;

MakeIndices = 1;
for L = 0:deltaL:evenLmax
  Lparity = parity(L);
  evenL = (Lparity==1);
  if (~evenL) & (L>oddLmax), continue; end    
  for jK = jkmn:2:1
    Kmx = min(L,Kmax);
    for K = 0:deltaK:Kmx
      evenK = parity(K)==1;
      if ((K==0) & (Lparity~=jK)), continue; end
      
      Mmx = min(L,Mmax);
      for M = -Mmx:Mmx


        for pS = pSmin:1
          qSmx = 1 - abs(pS);
          for qS = -qSmx:2:qSmx
            for pI = -pImax:pImax
              if ((~DirTilt)&((pI+pS-1)~=M)), continue; end
              
              qImax = 2*I - abs(pI);
              for qI = -qImax:2:qImax
                
                iRow = iRow + 1;
                if MakeIndices
                  L_(iRow) = L;
                  jK_(iRow) = jK;
                  K_(iRow) = K;
                  M_(iRow) = M;
                end
                
              end % qI
            end % pI
          end % qS
        end % pS

      end % M
    end % K
  end % jK
end % L

nRows = iRow;

Indices.L = L_;
Indices.jK = jK_;
Indices.K = K_;
Indices.M = M_;

return
%==================================================================

function p = parity(a)
if (mod(a,2)==0), p = +1; else p = -1; end
return
