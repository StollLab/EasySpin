% contfracspec   Compute spectrum using continued-fraction method
%
%   spec = contfracspec(z,alpha,beta,k)
%
%   z     ... abscissa of spectrum
%   alpha ... diagonal of tridiagonal matrix
%   beta  ... superdiagonal of tridiagonal matrix
%   k     ... number of elements to use

function spec = chili_contfracspec(z,alpha,beta,k,Method)

if (nargin<4), k = numel(alpha); end
if ((numel(alpha)<k) || (numel(beta)<k))
  error('alpha and/or beta have less than than %d elements!',k);
end

% Evaluate continued-fraction expression
% (Schneider/Freed 1989, p.34, Eq. 124)

if (nargin<5), Method = 0; end

switch Method
case 0
  % Right-to-left evaluation
  spec = inf;
  for m = k:-1:1
    spec = z + alpha(m) - beta(m)^2./spec;
  end
  spec = 1./spec;

case 1
  % Left-to-right evaluation:
  % Modified Lentz method (see Numerical Recipes)
  tiny = 1e-30;
  spec = tiny;
  C = spec;
  D = 0;
  for q = 1:k
    Cold = C;
    Dold = D;
    specold = spec;
    a = z + alpha(q);
    if (q==1), b = 1; else b = -beta(q-1)^2; end
    C = a + b./Cold;
    %C(C==0) = tiny;
    D = a + b.*Dold;
    %D(D==0) = tiny;
    D = 1./D;
    delta = C.*D;
    spec = specold.*delta;
  end
end

return
