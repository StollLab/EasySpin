% wignerd  Wigner D-matrix / D-function
%
%   D = wignerd(J,alpha,beta,gamma);
%   D = wignerd(J,alpha,beta,gamma,'-');
%   D = wignerd(J,alpha,beta,gamma,'+');
%   D = wignerd(J,beta,...);
%   Dm1m2 = wignerd(Jm1m2,...);
%
%   This function computes the Wigner rotation matrix with elements D^J_(m1,m2),
%   where m1 and m2 run from J to -J. If m1 and m2 are given in Jm1m2 = [J m1 m2],
%   only the single matrix element D^J_(m1,m2) is calculated.
%
%   J      ... rank J; J = 0, 1/2, 1, 3/2, 2, 5/2, etc.
%   Jm1m2  ... rank and two indices, [J m1 m2]; m1|m2 = J, J-1, ..., -J
%   alpha, beta, gamma ... Euler angles, in radians
%   D      ... Wigner rotation matrix, size (2J+1)x(2J+1)
%   Dm1m2  ... Wigner rotation matrix element
%
%   If only beta is given, alpha and gamma are assumed zero, and the Wigner
%   small-d matrix d^J(beta) or small-d matrix element d^J_(m1,m2)(beta)
%   is computed.
%
%   The sign convention for the rotation operator is given by the last argument
%   - '-': expm(-1i*alpha*Jz)*expm(-1i*beta*Jy)*expm(-1i*gamma*Jz)
%     (sign convention as in Brink/Satcher, Zare, Sakurai, Varshalovich).
%   - '+': expm(+1i*alpha*Jz)*expm(+1i*beta*Jy)*expm(+1i*gamma*Jz)
%     (sign convention as in Edmonds, Mathematica).
%   The default sign convention is '-'.
%
%   The basis is ordered +J..-J from left to right and from top to bottom,
%   so that e.g. the output matrix element D(1,2) corresponds to D^J_(J,J-1).

function D = wignerd(J,varargin)

defaultphase = '-';

switch nargin
  case 0
    help(mfilename);
    return
  case 1 % J
    error('Angles are missing');
  case 2 % J/Jm1m2, beta
    beta = varargin{1};
    betaonly = true;
    phase = defaultphase;
  case 3 % J/Jm1m2, beta, phase
    beta = varargin{1};
    betaonly = true;
    phase = varargin{2};
  case 4 % J/Jm1m2, alpha, beta, gamma
    alpha = varargin{1};
    beta = varargin{2};
    gamma = varargin{3};
    phase = defaultphase;
    betaonly = false;
  case 5 % J/Jm1m2, alpha, beta, gamma, phase
    alpha = varargin{1};
    beta = varargin{2};
    gamma = varargin{3};
    betaonly = false;
    phase = varargin{4};
end

switch numel(J)
  case 1
    calculateMatrix = true;
  case 3
    calculateMatrix = false;
    m1 = J(2);
    m2 = J(3);
    J = J(1);
  otherwise
    error('First input must be either J or [J m1 m2]');
end

if numel(J)~=1 || ~isreal(J) || mod(J,0.5) || J<0 || ~isnumeric(J)
  error('J must be one of 0, 1/2, 1, 3/2, 2, etc.');
end

if ~calculateMatrix
  if numel(m1)~=1 || ~isreal(m1) || m1<-J || m1>J
    error('m1 must be an integer between -J and J. Your J is %d.',J);
  end
  if numel(m2)~=1 || ~isreal(m2) || m2<-J || m2>J
    error('m1 must be an integer between -J and J. Your J is %d.',J);
  end
end

if ~ischar(phase) || numel(phase)~=1
  error('Last argument must be either ''+'' or ''-''.');
end

if (phase=='+')
  % do nothing
elseif (phase=='-')
  beta = -beta;
  if ~betaonly
    alpha = -alpha;
    gamma = -gamma;
  end
else
  error('Phase (last argument) must be either ''+'' or ''-''.');
end

if calculateMatrix
  % Calculate full Wigner D matrix
  %-----------------------------------------------------------------------------
  if ~any(beta) && (betaonly || (~any(alpha) && ~any(gamma)))
    D = eye(2*J+1);
    return;
  end
  
  % Calculate Wigner little-d matrix via Jy operator matrix exponential
  v = sqrt((1:2*J).*(2*J:-1:1))/2; % off-diagonals of Jy matrix (without 1/i)
  Jy = diag(v,+1) - diag(v,-1); % assemble Jy matrix
  d = expm(beta*Jy); % (i dropped since 1/i dropped in Jy)
  
  % Include alpha and gamma terms if given
  if betaonly
    D = d;
  else
    mz = J:-1:-J; % diagonal of Jz matrix
    D = (exp(1i*alpha*mz).'*exp(1i*gamma*mz)).*d;
  end
  
else
  % Calculate single matrix element of Wigner D matrix
  %-----------------------------------------------------------------------------
  % Calculate Wigner d-function via Jacobi polynomial
  % (see Biedenharn/Louck, Angular Momentum in Quantum Physics)
  
  % Determine k, a (nonnegative), and b (nonnegative)
  [k,idx] = min([J+m2,J-m2,J+m1,J-m1]);
  switch idx
    case 1, a = m1-m2; lam = 0;
    case 2, a = m2-m1; lam = a;
    case 3, a = m2-m1; lam = a;
    case 4, a = m1-m2; lam = 0;
  end
  b = 2*J-2*k-a;
  
  % Calculate quotient of binomial factors, nchoosek(2*J-k,k+a)/nchoose(k+b,k)
  q_enum = [2*J-k:-1:2*J-2*k-a+1, 1:b];
  q_denom = [1:k+a, k+1:k+b];
  if J > 28
    % Sort values in order not to loose accuracy
    q_enum = sort(q_enum);
    q_denom = sort(q_denom);
  end
  q = prod(q_enum./q_denom);
  
  % Calculate d-function
  d = (-1)^lam * sqrt(q) * ...
    sin(beta/2).^a .* cos(beta/2).^b .* ...
    jacobip(k,a,b,cos(beta));
  
  % Include alpha and gamma factors if given
  if betaonly
    D = d;
  else
    D = exp(1i*(alpha*m1+gamma*m2)).*d;
  end
  
end


%===============================================================================
% jacobip    Jacobi polynomial
%
%  P = jacobip(n,a,b,x)
%
% Calculate Jacobi polynomial of degree n with parameters a and b at x,
% utilizing the three-term recurrence relation with respect to n.
%-------------------------------------------------------------------------------

function P = jacobip(n,a,b,x)

P0 = 1;
if n==0, P = P0; return; end

P1 = (a-b)/2 + (1+(a+b)/2)*x;
if n==1, P = P1; return; end

a2b2 = a^2-b^2;
ab = a+b;
for N = 2:n
  P = ((2*N+ab-1)*((2*N+ab)*(2*N+ab-2)*x + a2b2).*P1 ...
    - 2*(N+a-1)*(N+b-1)*(2*N+ab)*P0)/(2*N*(N+ab)*(2*N+ab-2));
  P0 = P1;
  P1 = P;
end
