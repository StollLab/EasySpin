% plegendre    Legendre polynomials and Associated Legendre polynomials
%
%    y = plegendre(L,x)
%    y = plegendre(L,M,x)
%
%    Computes the Legendre polynomial of degree L, or the
%    associated Legendre polynomial of degree L and order M,
%    evaluated at x.
%
%    L and M are integers with L >= 0  and |M| <= L. x can
%    be a scalar or an array of real numbers with |x|<=1.

% For a faster and numerically better algorithm see
% Holmes, Featherstone,
% A unified approach to the Clenshaw summation and the recursive
% computation of very high degree and order normalised associated
% Legendre functions.
% Journal of Geodesy 2002 76:279-299
% https://doi.org/10.1007/s00190-002-0216-2

function y = plegendre(L,M,x)

if (nargin==0), help(mfilename); return; end

if (nargin==2)
  x = M;
  M = 0;
end

if (L<0) || (abs(M)>L) || any(abs(x(:))>1)
  error('Wrong Associated Legendre Polynomial parameters!');
end

%----------------------------------------------------------------
% Identities used in the computation
%  (see Press et al, Numerical Recipes in C++)
%----------------------------------------------------------------
%   M
%  P (z) == (2M-1)!! (1-z^2)^(M/2)            for L == M, M>=0
%   L
%
%      (-1)^M is the famous Condon-Shortley phase.
%
%   M                 M
%  P (z) == z (2M+1) P (x)                           for L == M+1, M>=0
%   L                 M
%
%
%   M        2N-1     M         N+M-1   M
%  P (z) == ------ z P (z)  -  ------- P (z)         for L >= M+2, M>=0
%   N        N-M      N-1        N-M    N-2
%
%----------------------------------------------------------------
%
%   -M              (N-M)!   M
%  P (z) == (-1)^M  ------  P (z)
%   N               (N+M)!   N
%----------------------------------------------------------------

negativeM = (M<0);
if (negativeM), M = -M; end

if (M>0)
  P0 = prod(1:2:2*M-1)*sqrt(1-x.^2).^M;
else
  P0 = ones(size(x));
end

if (L==M)
  y = P0;
else
  P1 = (2*M+1)*x.*P0;
  if (L==M+1)
    y = P1;
  else
    for n = M+2:L
      P2 = (2*n-1)/(n-M)*x.*P1 - (n+M-1)/(n-M)*P0;
      P0 = P1;
      P1 = P2;
    end
    y = P2;
  end
end

if (negativeM)
  y = (-1)^M/prod(L-M+1:L+M) * y;
end

return
