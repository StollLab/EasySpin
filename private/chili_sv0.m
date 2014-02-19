function v = chili_startingvector(Sys,Basis,Diffusion,Options)

if nargin<4, Options = struct('ununsed',NaN); end

if ~isfield(Options,'Tolerance')
  Options.Tolerance = 1e-10; % for numerical integration
end

if ~isfield(Options,'SeparateTransitions')
  Options.SeparateTransitions = 0;
end

evenLmax = Basis.evenLmax;
oddLmax = Basis.oddLmax;
Kmax = Basis.Kmax;
Mmax = Basis.Mmax;
deltaL = Basis.deltaL;
deltaK = Basis.deltaK;
jkmn = Basis.jKmin;
pSmin = Basis.pSmin;

DirTilt = (Sys.psi~=0);

Potential = any(Diffusion.lambda);
if (Potential)
  %AxialPotential = all(Diffusion.lambda([2 4])==0);
  PotLmax = size(Diffusion.xlk,1)-1;
  PotKmax = PotLmax;
  lambda = Diffusion.lambda;
end

if isfield(Basis,'Size')
  maxVals = Basis.Size;
else
  maxVals = 10000;
end

Value = zeros(maxVals,1);
Trans = zeros(maxVals,1);
Row   = zeros(maxVals,1);

nRow = 0; % number of rows of the column vector v
idx = 0;

iseven = @(x)mod(x,2)==0;
isodd = @(x)mod(x,2)~=0;
for L = 0:deltaL:evenLmax
  if isodd(L) && (L>oddLmax), continue; end    
  for jK = jkmn:2:1
    Kmx = min(L,Kmax);
    for K = 0:deltaK:Kmx
      if (K==0) && (parity(L)~=jK), continue; end
      
      %==============================================
      thisValue = 0;
      if iseven(L) && iseven(K)
        if (~Potential)
          % Zero potential: non-zero values for L==K==0
          if (L==0) && (K==0), thisValue = 1; end
        elseif ((PotKmax==0) && (K~=0))
          % Axial potential (max K == 0): zero for K>0
          %thisValue = 0;
        else
          %
          % Axial potential K==0, nonaxial potential: numerical integration
          thisValue = quadl(@orifun,0,1,Options.Tolerance,0,L,K,lambda);
          if (K~=0)
            thisValue = thisValue * sqrt((2*L+1)/prod(L-K+1:L+K));
          else
            thisValue = thisValue * sqrt((2*L+1)/2);
          end
        end
      end
      %==============================================
      
      Mmx = min(L,Mmax);
      for M = -Mmx:Mmx
        for pS = pSmin:1
          qSmx = 1 - abs(pS);
          for qS = -qSmx:2:qSmx
              if ((~DirTilt)&&((0+pS-1)~=M)), continue; end % Meirovich Eq.(A47)
              
              %==============================================
              NonZeroElement = (jK==1) && (M==0) && (pS~=0);
              %==============================================
                              
                nRow = nRow + 1;
                
                %==============================================
                if NonZeroElement && ...
                    (((~DirTilt) && (pS==1)) || ...
                    (( DirTilt) && (abs(pS)==1)))
                  if (thisValue~=0)
                    idx = idx + 1;
                    Value(idx) = thisValue;
                    Trans(idx) = 0; % gives mI of transition
                    Row(idx) = nRow;
                  end
                end
                %==============================================

          end
        end
      end
    end
  end
end

Value = Value(1:idx);
Trans = Trans(1:idx);
Row = Row(1:idx);

if (Options.SeparateTransitions)
  I = 0;
  v = full(sparse(Row,Trans+I+1,Value,nRow,2*I+1));
else
  v = full(sparse(Row,1,Value,nRow,1));
end

for iCol = 1:size(v,2)
  v(:,iCol) = v(:,iCol)/norm(v(:,iCol));
end

return
%==================================================================

function p = parity(a)
if mod(a,2), p = -1; else p = +1; end
return

%==================================================================
% orifun  Integrand of the orientational integral in the starting vector.
% see Schneider 1989, p.19
function val = orifun(z,L,K,lambda)

% based on the integral (valid for even K)
% \int_0^{2\pi} cos(K g) exp(B cos(2 g)) d g =
% 2 \pi besseli(K/2,B)

% Potential terms with K == 0
A = 0;
if lambda(1), A = A + lambda(1)/2*plegendre(2,0,z); end
if lambda(3), A = A + lambda(3)/2*plegendre(4,0,z); end

% Potential terms with K == 2
B = 0;
if lambda(2), B = B + lambda(2)/(2*sqrt( 6))*plegendre(2,2,z); end
if lambda(4), B = B + lambda(4)/(6*sqrt(10))*plegendre(4,2,z); end

val = plegendre(L,K,z).*exp(A).*besseli(K/2,B)*(2*pi);

return
%====================================================================
