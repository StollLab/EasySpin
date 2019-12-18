% contfracspec   Compute spectrum using continued-fraction method
%
%   spec = contfracspec(z,alpha,beta,k)
%   spec = contfracspec(z,alpha,beta,k,useLentz)
%
% Inputs:
%   z        ... abscissa of spectrum
%   alpha    ... diagonal of tridiagonal matrix
%   beta     ... superdiagonal of tridiagonal matrix
%   k        ... number of elements to use
%   useLentz ... if true, use modified Lentz method; default is false
%
% Outputs:
%   spe   c  ... calculated spectrum

% Evaluate continued-fraction expression
% (seee Schneider/Freed, Biol. Magn. Reson. 1989, p.34, Eq. 124)

function spec = chili_contfracspec(z,alpha,beta,k,useLentz)

% Supplement optional inputs, input validation
if nargin<4, k = numel(alpha); end
if numel(alpha)<k || numel(beta)<k
  error('alpha and/or beta have less than than %d elements!',k);
end
if nargin<5, useLentz = false; end

if useLentz
  
  % Left-to-right evaluation:
  % Modified Lentz method (see Press et al, Numerical Recipes in C, 2nd ed., section 5.2)
  tiny = 1e-30;
  spec = tiny;
  C = spec;
  D = 0;
  for q = 1:k
    Cold = C;
    Dold = D;
    specold = spec;
    b = z + alpha(q);
    if q==1, a = 1; else, a = -beta(q-1)^2; end
    C = b + a./Cold;
    %C(C==0) = tiny;
    D = b + a.*Dold;
    %D(D==0) = tiny;
    D = 1./D;
    Delta = C.*D;
    spec = specold.*Delta;
  end
  
else
  
  % Right-to-left evaluation, straightforward
  spec = inf;
  for m = k:-1:1
    spec = z + alpha(m) - beta(m)^2./spec;
  end
  spec = 1./spec;
  
end

return
