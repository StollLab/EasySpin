% Build basis with basis functions and ordering as used in the Freed program
% (S = 1/2 and zero, one, or two nuclei).
function Basis = chili_basisbuild(Basis,Sys)

DirTilt = Basis.DirTilt;

evenLmax = Basis.evenLmax;
oddLmax = Basis.oddLmax;

jKmin = Basis.jKmin;
Kmax = Basis.Kmax;
deltaK = Basis.evenK+1;

Mmax = Basis.Mmax;

pSmin = Basis.pSmin;
pImax = Basis.pImax;
pI1max = pImax; pI2max = pImax;

MpSymm = Basis.MpSymm;

I = Sys.I;
nNuclei = numel(I);
if nNuclei>=1, I1 = I(1); end
if nNuclei>=2, I2 = I(2); end

iRow = 0;
iSpatial = 0;

makeIndices = true;

if makeIndices
  nRowBlock = 20000; % size of pre-allocation block
  Indices = zeros(nRowBlock,4+2+2*nNuclei);
end

parity = @(x) 1 - 2*mod(x,2); % -1 for odd number, +1 for even number

for L = 0:evenLmax
  Lparity = parity(L);
  evenL = (Lparity==1);
  if (~evenL) && (L>oddLmax), continue; end    
  for jK = jKmin:2:1
    for K = 0:deltaK:min(L,Kmax)
      if (K==0) && (Lparity~=jK), continue; end
      Mmx = min(L,Mmax);
      for M = -Mmx:Mmx
        iSpatial = iSpatial + 1;
        
        
        for pS = pSmin:1
          qSmx = 1 - abs(pS);
          for qS = -qSmx:2:qSmx
            
            % no nuclei ---------------------------------------
            if nNuclei==0
              if MpSymm && ~DirTilt && (pS-1)~=M, continue; end % Meirovitch Eq.(A47)
              
              iRow = iRow + 1;
              if makeIndices
                if mod(iRow,nRowBlock)==0
                  Indices(iRow+nRowBlock,:) = 0;
                end
                Indices(iRow,:) = [L jK K M pS qS];
              end
                
            % one nucleus ---------------------------------------
            elseif nNuclei==1
              
              for pI1 = -pI1max:pI1max
                if MpSymm && ~DirTilt && pI1+pS-1~=M, continue; end % Meirovich Eq.(A47)
                qI1max = 2*I1 - abs(pI1);
                for qI1 = -qI1max:2:qI1max
                  
                  iRow = iRow + 1;
                  if makeIndices
                    if mod(iRow,nRowBlock)==0
                      Indices(iRow+nRowBlock,:) = 0;
                    end
                    Indices(iRow,:) = [L jK K M pS qS pI1 qI1];
                  end
                  
                end % qI
              end % pI
              
            % two nuclei ---------------------------------------
            elseif nNuclei==2
              
              for pI1 = -pI1max:pI1max
                qI1max = 2*I1 - abs(pI1);
                for qI1 = -qI1max:2:qI1max
                  for pI2 = -pI2max:pI2max
                    if MpSymm && ~DirTilt && pI1+pI2+pS-1~=M, continue; end % Meirovich Eq.(A47)
                    qI2max = 2*I2 - abs(pI2);
                    for qI2 = -qI2max:2:qI2max
                      
                      iRow = iRow + 1;
                      if makeIndices
                        if mod(iRow,nRowBlock)==0
                          Indices(iRow+nRowBlock,:) = 0;
                        end
                        Indices(iRow,:) = [L jK K M pS qS pI1 qI1 pI2 qI2];
                      end
                      
                    end % qI2
                  end % pI2
                end % qI1
              end % pI1

            end
            
          end % qS
        end % pS

      end % M
    end % K
  end % jK
end % L
Indices  = Indices(1:iRow,:);

Basis.L = Indices(:,1);
Basis.jK = Indices(:,2);
Basis.K = Indices(:,3);
Basis.M = Indices(:,4);
Basis.pS = Indices(:,5);
Basis.qS = Indices(:,6);
if nNuclei>=1
  Basis.pI1 = Indices(:,7);
  Basis.qI1 = Indices(:,8);
  if nNuclei>=2
    Basis.pI2 = Indices(:,9);
    Basis.qI2 = Indices(:,10);
  end
end

return
