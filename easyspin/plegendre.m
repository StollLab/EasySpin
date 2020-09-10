% plegendre    Legendre polynomials and associated Legendre polynomials
%
%    y = plegendre(L,z)
%    y = plegendre(L,M,z)
%    y = plegendre(L,M,z,CSphase)
%
%    Computes the Legendre polynomial of degree L, or the associated Legendre
%    polynomial of degree L and order M, evaluated at z.
%
%    L and M are integers with L >= 0  and |M| <= L.
%    z can be a scalar or an array of real numbers with -1 <= z <= 1.
%
%    The optional input CSphase specifies whether the Condon-Shortley phase
%    (-1)^M should be included (true) or not (false). The default is true.

% For a faster and numerically better algorithm see
%   S.A.Holmes, W.E.Featherstone,
%   A unified approach to the Clenshaw summation and the recursive
%   computation of very high degree and order normalised associated
%   Legendre functions.
%   Journal of Geodesy 2002 76:279-299
%   https://doi.org/10.1007/s00190-002-0216-2

%-------------------------------------------------------------------------------
% Identities used in the computation
%  (see Press et al, Numerical Recipes in C++)
%  (without the Condon-Shortley phase)
%-------------------------------------------------------------------------------
%   M
%  P (z) == (2M-1)!! (1-z^2)^(M/2)                   for L=M, M>=0
%   L
%
%   M                 M
%  P (z) == z (2M+1) P (x)                           for L=M+1, M>=0
%   L                 M
%
%
%   M        2L-1     M         L+M-1   M
%  P (z) == ------ z P (z)  -  ------- P (z)         for L>=M+2, M>=0
%   L        L-M      L-1        L-M    L-2
%
%-------------------------------------------------------------------------------
%
%   -M              (L-M)!   M
%  P (z) == (-1)^M  ------  P (z)
%   L               (L+M)!   L
%
%-------------------------------------------------------------------------------

function y = plegendre(varargin)

if nargin==0, help(mfilename); return; end

if nargin<2
  error('At least 2 inputs are required (L, z).');
end

if nargin>4
  error('At most 4 inputs are possible (L, M, z, CSphase).');
end

if nargin<4
  CondonShortleyPhase = true;
end

switch nargin
  case 2
    L = varargin{1};
    M = 0;
    z = varargin{2};
  case 3
    L = varargin{1};
    M = varargin{2};
    z = varargin{3};
  case 4
    L = varargin{1};
    M = varargin{2};
    z = varargin{3};
    CondonShortleyPhase = varargin{4};
end

if numel(L)~=1 || L<0 || mod(L,1)
  error('L (first input) must be a non-negative integer (0, 1, 2, etc).');
end
if numel(M)~=1 || abs(M)>L || mod(M,1)
  error('M (second input) must be an integer with -L<=M<=L.');
end
if ~isreal(z) || any(abs(z(:))>1)
  error('All values in z (third input) must satisfy -1<=z<=1.');
end
if numel(CondonShortleyPhase)~=1 || ~islogical(CondonShortleyPhase)
  error('CSphase (4th input) must be true or false.');
end

% Calculate P_L^M for non-negative M first.
negativeM = M<0;
if negativeM
  M = -M;
end

% Calculate P_M^M
if M>0
  P0 = prod(1:2:2*M-1)*sqrt(1-z.^2).^M;
else
  P0 = ones(size(z));
end

% Apply recurrcence relation to increase degree from M to L
if L==M
  y = P0;
else
  P1 = (2*M+1)*z.*P0;
  if L==M+1
    y = P1;
  else
    for L_ = M+2:L
      P2 = (2*L_-1)/(L_-M)*z.*P1 - (L_+M-1)/(L_-M)*P0;
      P0 = P1;
      P1 = P2;
    end
    y = P2;
  end
end

% If M was negative, compute P for negative M
if negativeM
  % (The (-1)^M in this expression is not the Condon-Shortley phase.)
  y = (-1)^M/prod(L-M+1:L+M) * y;
end

% Include Condon-Shortley phase if requested
if CondonShortleyPhase
  y = (-1)^M * y;
end

return
