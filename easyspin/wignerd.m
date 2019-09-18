% wignerd  Wigner rotation matrices
%
%   D = wignerd(J,alpha,beta,gamma);
%   D = wignerd(J,alpha,beta,gamma,phase);
%   Dm1m2 = wignerd(Jm1m2,alpha,beta,gamma);
%   Dm1m2 = wignerd(Jm1m2,alpha,beta,gamma,phase);
%   d = wignerd(J,beta);
%   d = wignerd(J,beta,phase);
%   dm1m2 = wignerd(Jm1m2,beta);
%   dm1m2 = wignerd(Jm1m2,beta,phase);
%
%   This function computes the Wigner rotation matrix with elements D^J_{m1,m2},
%   where m1 and m2 run from J to -J. If m1 and m2 are given in Jm1m2 = [J m1 m2],
%   only the single matrix element D^J_{m1,m2} is calculated.
%
%   Input:
%     J      ... rank J; J = 0, 1/2, 1, 3/2, 2, 5/2, etc.
%     Jm1m2  ... rank and two indices, [J m1 m2]; with m1,m2 = J, J-1, ..., -J
%     alpha, beta, gamma ... Euler angles, in radians
%             If only beta is given, alpha and gamma are assumed zero.
%     phase  ... sign convention for the rotation operator, '+' or '-'
%            '-': expm(-1i*alpha*Jz)*expm(-1i*beta*Jy)*expm(-1i*gamma*Jz)
%                (as used in Brink/Satcher, Zare, Sakurai, Varshalovich, Biedenharn/Louck, Mehring).
%            '+': expm(+1i*alpha*Jz)*expm(+1i*beta*Jy)*expm(+1i*gamma*Jz)
%                (as used in Edmonds, Mathematica).
%            The default sign convention is '-'.
%
%  Output: 
%     D      ... Wigner rotation matrix D^J, size (2J+1)x(2J+1)
%     Dm1m2  ... Wigner rotation matrix element D^J_{m1,m2}
%     d      ... reduced Wigner rotation matrix d^J, size (2J+1)x(2J+1)
%     dm1m2  ... reduced Wigner rotation matrix element d^J_{m1,m2}
%
%   The basis is ordered +J..-J from left to right and from top to bottom,
%   so that e.g. the output matrix element D(1,2) corresponds to D^J_{J,J-1}.

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

% Determine phase setting (convert from char to +-1)
if ~ischar(phase) || numel(phase)~=1
  error('Last argument must be either ''+'' or ''-''.');
end
if phase=='-'
  phase = -1;
else
  phase = +1;
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
  d = expm(phase*beta*Jy); % (i dropped since 1/i dropped in Jy)
  
  % Include alpha and gamma terms if given
  if betaonly
    D = d;
  else
    mz = J:-1:-J; % diagonal of Jz matrix
    D = (exp(1i*phase*alpha*mz).'*exp(1i*phase*gamma*mz)).*d;
  end
  
else
  
  % Calculate single matrix element of Wigner D matrix
  %-----------------------------------------------------------------------------
  
  % Calculate little-d function of beta angle
  if J<=4
    d = littled_explicit(J,m1,m2,phase*beta);
  else
    d = littled_jacobi(J,m1,m2,phase*beta);
  end
  
  % Include alpha and gamma factors if given
  if betaonly
    D = d;
  else
    D = exp(phase*1i*(alpha*m1+gamma*m2)).*d;
  end
  
end


%===============================================================================
% jacobip    Jacobi polynomial
%
%  P = jacobip(n,a,b,x)
%
% Calculate Jacobi polynomial of degree n with parameters a and b at x,
% utilizing the three-term recurrence relation with respect to n.
% see Eqs. 18.9.1 and 18.9.2 at https://dlmf.nist.gov/18.9
% see http://functions.wolfram.com/Polynomials/JacobiP/17/01/01/01/0002/
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

%===============================================================================
% Calculate reduced rotation matrix element d^J_{m1,m2} via Jacobi polynomial
% see
%    L.C.Biedenharn & J.D. Louck, Angular Momentum in Quantum Physics,
%      eq. 3.73, 3.74, p.50
%    Varshalovich et al, Quantum Theory of Angular Momentum
%      eq. 4.3.4 (13)-(15)
%-------------------------------------------------------------------------------
function d = littled_jacobi(J,m1,m2,beta)

% Invert sign, since the Biedenharn/Louck equations are for expm(-1i*beta*Jy)
beta = -beta;

% Determine k, a (nonnegative), and b (nonnegative)
[k,idx] = min([J+m2,J-m2,J+m1,J-m1]);
switch idx
  case 1, mu = m1-m2; nu = -m1-m2; lam = mu;
  case 2, mu = m2-m1; nu =  m1+m2; lam = 0;
  case 3, mu = m2-m1; nu = -m1-m2; lam = 0;
  case 4, mu = m1-m2; nu =  m1+m2; lam = mu;
end

% Calculate quotient of binomial factors, nchoosek(2*J-k,k+a)/nchoose(k+b,k)
q_enum = [2*J-k:-1:2*J-2*k-mu+1, 1:nu];
q_denom = [1:k+mu, k+1:k+nu];
if J > 28
  % Sort values in order not to loose accuracy
  q_enum = sort(q_enum);
  q_denom = sort(q_denom);
end
q = prod(q_enum./q_denom);

% Calculate d-function
d = (-1)^lam * sqrt(q) * ...
  sin(beta/2).^mu .* cos(beta/2).^nu .* jacobip(k,mu,nu,cos(beta));

%===============================================================================
% Calculate reduced rotation matrix element d^J_{m1,m2} using explicit expressions
% (see Varshalovich et al, Tables 4.3, 4.4, 4.5, and 4.6, p. 119)
%-------------------------------------------------------------------------------
function d = littled_explicit(J,m1,m2,beta)

% Invert sign, since the expressions from Varshalovich, used below, are derived
% for expm(-1i*Jy*beta).
beta = -beta;

ridx = J+1-m1;
cidx = J+1-m2;
linidx = (cidx-1)*(2*J+1)+ridx;

switch J
  case 0
    d = ones(size(beta));
  case 0.5
    switch linidx
      case 1, d = cos(beta/2);
      case 2, d = sin(beta/2);
      case 3, d = -sin(beta/2);
      case 4, d = cos(beta/2);
    end
  case 1
    switch linidx
      case 1, d = (1+cos(beta))/2;
      case 2, d = sin(beta)/sqrt(2);
      case 3, d = (1-cos(beta))/2;
      case 4, d = -sin(beta)/sqrt(2);
      case 5, d  = cos(beta);
      case 6, d = sin(beta)/sqrt(2);
      case 7, d = (1-cos(beta))/2;
      case 8, d = -sin(beta)/sqrt(2);
      case 9, d = (1+cos(beta))/2;
    end
  case 1.5
    switch linidx
      case  1, d = cos(beta/2).^3;
      case  2, d = sqrt(3)*sin(beta/2).*cos(beta/2).^2;
      case  3, d = sqrt(3)*sin(beta/2).^2.*cos(beta/2);
      case  4, d = sin(beta/2).^3;
      case  5, d = -sqrt(3)*sin(beta/2).*cos(beta/2).^2;
      case  6, d = cos(beta/2).*(3*cos(beta/2).^2-2);
      case  7, d = -sin(beta/2).*(3*sin(beta/2).^2-2);
      case  8, d = sqrt(3)*sin(beta/2).^2.*cos(beta/2);
      case  9, d = sqrt(3)*sin(beta/2).^2.*cos(beta/2);
      case 10, d = sin(beta/2).*(3*sin(beta/2).^2-2);
      case 11, d = cos(beta/2).*(3*cos(beta/2).^2-2);
      case 12, d = sqrt(3)*sin(beta/2).*cos(beta/2).^2;
      case 13, d = -sin(beta/2).^3;
      case 14, d = sqrt(3)*sin(beta/2).^2.*cos(beta/2);
      case 15, d = -sqrt(3)*sin(beta/2).*cos(beta/2).^2;
      case 16, d = cos(beta/2).^3;
    end
  case 2
    switch linidx
      case  1, d = (1+cos(beta)).^2/4;
      case  2, d = sin(beta).*(1+cos(beta))/2;
      case  3, d = sqrt(3/2)/2*sin(beta).^2;
      case  4, d = sin(beta).*(1-cos(beta))/2;
      case  5, d = (1-cos(beta)).^2/4;
      case  6, d = -sin(beta).*(1+cos(beta))/2;
      case  7, cb = cos(beta); d = (2*cb.^2+cb-1)/2;
      case  8, d = sqrt(3/2)*sin(beta).*cos(beta);
      case  9, cb = cos(beta); d = -(2*cb.^2-cb-1)/2;
      case 10, d = sin(beta).*(1-cos(beta))/2;
      case 11, d = sqrt(3/2)/2*sin(beta).^2;
      case 12, d = -sqrt(3/2)*sin(beta).*cos(beta);
      case 13, d = (3*cos(beta).^2-1)/2;
      case 14, d = sqrt(3/2)*sin(beta).*cos(beta);
      case 15, d = sqrt(3/2)/2*sin(beta).^2;
      case 16, d = -sin(beta).*(1-cos(beta))/2;
      case 17, cb = cos(beta); d = -(2*cb.^2-cb-1)/2;
      case 18, d = - sqrt(3/2)*sin(beta).*cos(beta);
      case 19, cb = cos(beta); d = (2*cb.^2+cb-1)/2;
      case 20, d = sin(beta).*(1+cos(beta))/2;
      case 21, d = (1-cos(beta)).^2/4;
      case 22, d = -sin(beta).*(1-cos(beta))/2;
      case 23, d = sqrt(3/2)/2*sin(beta).^2;
      case 24, d = -sin(beta).*(1+cos(beta))/2;
      case 25, d = (1+cos(beta)).^2/4;
    end
  case 2.5
    switch linidx
      case  1, d = cos(beta/2).^5;
      case  2, d = sqrt(5)*sin(beta/2).*cos(beta/2).^4;
      case  3, d = sqrt(10)*sin(beta/2).^2.*cos(beta/2).^3;
      case  4, d = sqrt(10)*sin(beta/2).^3.*cos(beta/2).^2;
      case  5, d = sqrt(5)*sin(beta/2).^4.*cos(beta/2);
      case  6, d = sin(beta/2).^5;
      case  7, d = -sqrt(5)*sin(beta/2).*cos(beta/2).^4;
      case  8, d = cos(beta/2).^3.*(1-5*sin(beta/2).^2);
      case  9, d = sqrt(2)*sin(beta/2).*cos(beta/2).^2.*(2-5*sin(beta/2).^2);
      case 10, d = -sqrt(2)*sin(beta/2).^2.*cos(beta/2).*(2-5*cos(beta/2).^2);
      case 11, d = -sin(beta/2).^3.*(1-5*cos(beta/2).^2);
      case 12, d = sqrt(5)*sin(beta/2).^4.*cos(beta/2);
      case 13, d = sqrt(10)*sin(beta/2).^2.*cos(beta/2).^3;
      case 14, d = -sqrt(2)*sin(beta/2).*cos(beta/2).^2.*(2-5*sin(beta/2).^2);
      case 15, d = cos(beta/2).*(3-12*cos(beta/2).^2+10*cos(beta/2).^4);
      case 16, d = sin(beta/2).*(3-12*sin(beta/2).^2+10*sin(beta/2).^4);
      case 17, d = -sqrt(2)*sin(beta/2).^2.*cos(beta/2).*(2-5*cos(beta/2).^2);
      case 18, d = sqrt(10)*sin(beta/2).^3.*cos(beta/2).^2;
      case 19, d = -sqrt(10)*sin(beta/2).^3.*cos(beta/2).^2;
      case 20, d = -sqrt(2)*sin(beta/2).^2.*cos(beta/2).*(2-5*cos(beta/2).^2);
      case 21, d = -sin(beta/2).*(3-12*sin(beta/2).^2+10*sin(beta/2).^4);
      case 22, d = cos(beta/2).*(3-12*cos(beta/2).^2+10*cos(beta/2).^4);
      case 23, d = sqrt(2)*sin(beta/2).*cos(beta/2).^2.*(2-5*sin(beta/2).^2);
      case 24, d = sqrt(10)*sin(beta/2).^2.*cos(beta/2).^3;
      case 25, d = sqrt(5)*sin(beta/2).^4.*cos(beta/2);
      case 26, d = sin(beta/2).^3.*(1-5*cos(beta/2).^2);
      case 27, d = -sqrt(2)*sin(beta/2).^2.*cos(beta/2).*(2-5*cos(beta/2).^2);
      case 28, d = -sqrt(2)*sin(beta/2).*cos(beta/2).^2.*(2-5*sin(beta/2).^2);
      case 29, d = cos(beta/2).^3.*(1-5*sin(beta/2).^2);
      case 30, d = sqrt(5)*sin(beta/2).*cos(beta/2).^4;
      case 31, d = -sin(beta/2).^5;
      case 32, d = sqrt(5)*sin(beta/2).^4.*cos(beta/2);
      case 33, d = -sqrt(10)*sin(beta/2).^3.*cos(beta/2).^2;
      case 34, d = sqrt(10)*sin(beta/2).^2.*cos(beta/2).^3;
      case 35, d = -sqrt(5)*sin(beta/2).*cos(beta/2).^4;
      case 36, d = cos(beta/2).^5;
    end
  case 3
    switch linidx
      case  1, d = (1+cos(beta)).^3/8;
      case  2, d = sqrt(6)/8*sin(beta).*(1+cos(beta)).^2;
      case  3, d = sqrt(15)/8*sin(beta).^2.*(1+cos(beta));
      case  4, d = sqrt(5)/4*sin(beta).^3;
      case  5, d = sqrt(15)/8*sin(beta).^2.*(1-cos(beta));
      case  6, d = sqrt(6)/8*sin(beta).*(1-cos(beta)).^2;
      case  7, d = (1-cos(beta)).^3/8;
      case  8, d = -sqrt(6)/8*sin(beta).*(1+cos(beta)).^2;
      case  9, d = -1/4*(1+cos(beta)).^2.*(2-3*cos(beta));
      case 10, d = -sqrt(10)/8*sin(beta).*(1-2*cos(beta)-3*cos(beta).^2);
      case 11, d = sqrt(30)/4*sin(beta).^2.*cos(beta);
      case 12, d = sqrt(10)/8*sin(beta).*(1+2*cos(beta)-3*cos(beta).^2);
      case 13, d = 1/4*(1-cos(beta)).^2.*(2+3*cos(beta));
      case 14, d = sqrt(6)/8*sin(beta).*(1-cos(beta)).^2;
      case 15, d = sqrt(15)/8*sin(beta).^2.*(1+cos(beta));
      case 16, d = sqrt(10)/8*sin(beta).*(1-2*cos(beta)-3*cos(beta).^2);
      case 17, d = -1/8*(1+cos(beta)).*(1+10*cos(beta)-15*cos(beta).^2);
      case 18, d = -sqrt(3)/4*sin(beta).*(1-5*cos(beta).^2);
      case 19, d = -1/8*(1-cos(beta)).*(1-10*cos(beta)-15*cos(beta).^2);
      case 20, d = sqrt(10)/8*sin(beta).*(1+2*cos(beta)-3*cos(beta).^2);
      case 21, d = sqrt(15)/8*sin(beta).^2.*(1-cos(beta));
      case 22, d = -sqrt(5)/4*sin(beta).^3;
      case 23, d = sqrt(30)/4*sin(beta).^2.*cos(beta);
      case 24, d = sqrt(3)/4*sin(beta).*(1-5*cos(beta).^2);
      case 25, d = -1/2*cos(beta).*(3-5*cos(beta).^2);
      case 26, d = -sqrt(3)/4*sin(beta).*(1-5*cos(beta).^2);
      case 27, d = sqrt(30)/4*sin(beta).^2.*cos(beta);
      case 28, d = sqrt(5)/4*sin(beta).^3;
      case 29, d = sqrt(15)/8*sin(beta).^2.*(1-cos(beta));
      case 30, d = -sqrt(10)/8*sin(beta).*(1+2*cos(beta)-3*cos(beta).^2);
      case 31, d = -1/8*(1-cos(beta)).*(1-10*cos(beta)-15*cos(beta).^2);
      case 32, d = sqrt(3)/4*sin(beta).*(1-5*cos(beta).^2);
      case 33, d = -1/8*(1+cos(beta)).*(1+10*cos(beta)-15*cos(beta).^2);
      case 34, d = -sqrt(10)/8*sin(beta).*(1-2*cos(beta)-3*cos(beta).^2);
      case 35, d = sqrt(15)/8*sin(beta).^2.*(1+cos(beta));
      case 36, d = -sqrt(6)/8*sin(beta).*(1-cos(beta)).^2;
      case 37, d = 1/4*(1-cos(beta)).^2.*(2+3*cos(beta));
      case 38, d = -sqrt(10)/8*sin(beta).*(1+2*cos(beta)-3*cos(beta).^2);
      case 39, d = sqrt(30)/4*sin(beta).^2.*cos(beta);
      case 40, d = sqrt(10)/8*sin(beta).*(1-2*cos(beta)-3*cos(beta).^2);
      case 41, d = -1/4*(1+cos(beta)).^2.*(2-3*cos(beta));
      case 42, d = sqrt(6)/8*sin(beta).*(1+cos(beta)).^2;
      case 43, d = (1-cos(beta)).^3/8;
      case 44, d = -sqrt(6)/8*sin(beta).*(1-cos(beta)).^2;
      case 45, d = sqrt(15)/8*sin(beta).^2.*(1-cos(beta));
      case 46, d = -sqrt(5)/4*sin(beta).^3;
      case 47, d = sqrt(15)/8*sin(beta).^2.*(1+cos(beta));
      case 48, d = -sqrt(6)/8*sin(beta).*(1+cos(beta)).^2;
      case 49, d = (1+cos(beta)).^3/8;
    end
  case 3.5
    switch linidx
      case  1, d = cos(beta/2).^7;
      case  2, d = sqrt(7)*cos(beta/2).^6.*sin(beta/2);
      case  3, d = sqrt(21)*cos(beta/2).^5.*sin(beta/2).^2;
      case  4, d = sqrt(35)*cos(beta/2).^4.*sin(beta/2).^3;
      case  5, d = sqrt(35)*cos(beta/2).^3.*sin(beta/2).^4;
      case  6, d = sqrt(21)*cos(beta/2).^2.*sin(beta/2).^5;
      case  7, d = sqrt(7)*cos(beta/2).*sin(beta/2).^6;
      case  8, d = sin(beta/2).^7;
      case  9, d = -sqrt(7)*cos(beta/2).^6.*sin(beta/2);
      case 10, d = cos(beta/2).^5.*(1-7*sin(beta/2).^2);
      case 11, d = sqrt(3)*cos(beta/2).^4.*sin(beta/2).*(2-7*sin(beta/2).^2);
      case 12, d = sqrt(5)*cos(beta/2).^3.*sin(beta/2).^2.*(3-7*sin(beta/2).^2);
      case 13, d = -sqrt(5)*cos(beta/2).^2.*sin(beta/2).^3.*(3-7*cos(beta/2).^2);
      case 14, d = -sqrt(3)*cos(beta/2).*sin(beta/2).^4.*(2-7*cos(beta/2).^2);
      case 15, d = -sin(beta/2).^5.*(1-7*cos(beta/2).^2);
      case 16, d = sqrt(7)*cos(beta/2).*sin(beta/2).^6;
      case 17, d = sqrt(21)*cos(beta/2).^5.*sin(beta/2).^2;
      case 18, d = -sqrt(3)*cos(beta/2).^4.*sin(beta/2).*(2-7*sin(beta/2).^2);
      case 19, d = cos(beta/2).^3.*(10-30*cos(beta/2).^2+21*cos(beta/2).^4);
      case 20, d = sqrt(15)*cos(beta/2).^2.*sin(beta/2).*(2-8*cos(beta/2).^2+7*cos(beta/2).^4);
      case 21, d = sqrt(15)*cos(beta/2).*sin(beta/2).^2.*(2-8*sin(beta/2).^2+7*sin(beta/2).^4);
      case 22, d = sin(beta/2).^3.*(10-30*sin(beta/2).^2+21*sin(beta/2).^4);
      case 23, d = -sqrt(3)*cos(beta/2).*sin(beta/2).^4.*(2-7*cos(beta/2).^2);
      case 24, d = sqrt(21)*cos(beta/2).^2.*sin(beta/2).^5;
      case 25, d = -sqrt(35)*cos(beta/2).^4.*sin(beta/2).^3;
      case 26, d = sqrt(5)*cos(beta/2).^3.*sin(beta/2).^2.*(3-7*sin(beta/2).^2);
      case 27, d = -sqrt(15)*cos(beta/2).^2.*sin(beta/2).*(2-8*cos(beta/2).^2+7*cos(beta/2).^4);
      case 28, d = -cos(beta/2).*(4-30*cos(beta/2).^2+60*cos(beta/2).^4-35*cos(beta/2).^6);
      case 29, d = sin(beta/2).*(4-30*sin(beta/2).^2+60*sin(beta/2).^4-35*sin(beta/2).^6);
      case 30, d = sqrt(15)*cos(beta/2).*sin(beta/2).^2.*(2-8*sin(beta/2).^2+7*sin(beta/2).^4);
      case 31, d = -sqrt(5)*cos(beta/2).^2.*sin(beta/2).^3.*(3-7*cos(beta/2).^2);
      case 32, d = sqrt(35)*cos(beta/2).^3.*sin(beta/2).^4;
      case 33, d = sqrt(35)*cos(beta/2).^3.*sin(beta/2).^4;
      case 34, d = sqrt(5)*cos(beta/2).^2.*sin(beta/2).^3.*(3-7*cos(beta/2).^2);
      case 35, d = sqrt(15)*cos(beta/2).*sin(beta/2).^2.*(2-8*sin(beta/2).^2+7*sin(beta/2).^4);
      case 36, d = -sin(beta/2).*(4-30*sin(beta/2).^2+60*sin(beta/2).^4-35*sin(beta/2).^6);
      case 37, d = -cos(beta/2).*(4-30*cos(beta/2).^2+60*cos(beta/2).^4-35*cos(beta/2).^6);
      case 38, d = sqrt(15)*cos(beta/2).^2.*sin(beta/2).*(2-8*cos(beta/2).^2+7*cos(beta/2).^4);
      case 39, d = sqrt(5)*cos(beta/2).^3.*sin(beta/2).^2.*(3-7*sin(beta/2).^2);
      case 40, d = sqrt(35)*cos(beta/2).^4.*sin(beta/2).^3;
      case 41, d = -sqrt(21)*cos(beta/2).^2.*sin(beta/2).^5;
      case 42, d = -sqrt(3)*cos(beta/2).*sin(beta/2).^4.*(2-7*cos(beta/2).^2);
      case 43, d = -sin(beta/2).^3.*(10-30*sin(beta/2).^2+21*sin(beta/2).^4);
      case 44, d = sqrt(15)*cos(beta/2).*sin(beta/2).^2.*(2-8*sin(beta/2).^2+7*sin(beta/2).^4);
      case 45, d = -sqrt(15)*cos(beta/2).^2.*sin(beta/2).*(2-8*cos(beta/2).^2+7*cos(beta/2).^4);
      case 46, d = cos(beta/2).^3.*(10-30*cos(beta/2).^2+21*cos(beta/2).^4);
      case 47, d = sqrt(3)*cos(beta/2).^4.*sin(beta/2).*(2-7*sin(beta/2).^2);
      case 48, d = sqrt(21)*cos(beta/2).^5.*sin(beta/2).^2;
      case 49, d = sqrt(7)*cos(beta/2).*sin(beta/2).^6;
      case 50, d = sin(beta/2).^5.*(1-7*cos(beta/2).^2);
      case 51, d = -sqrt(3)*cos(beta/2).*sin(beta/2).^4.*(2-7*cos(beta/2).^2);
      case 52, d = sqrt(5)*cos(beta/2).^2.*sin(beta/2).^3.*(3-7*cos(beta/2).^2);
      case 53, d = sqrt(5)*cos(beta/2).^3.*sin(beta/2).^2.*(3-7*sin(beta/2).^2);
      case 54, d = -sqrt(3)*cos(beta/2).^4.*sin(beta/2).*(2-7*sin(beta/2).^2);
      case 55, d = cos(beta/2).^5.*(1-7*sin(beta/2).^2);
      case 56, d = sqrt(7)*cos(beta/2).^6.*sin(beta/2);
      case 57, d = -sin(beta/2).^7;
      case 58, d = sqrt(7)*cos(beta/2).*sin(beta/2).^6;
      case 59, d = -sqrt(21)*cos(beta/2).^2.*sin(beta/2).^5;
      case 60, d = sqrt(35)*cos(beta/2).^3.*sin(beta/2).^4;
      case 61, d = -sqrt(35)*cos(beta/2).^4.*sin(beta/2).^3;
      case 62, d = sqrt(21)*cos(beta/2).^5.*sin(beta/2).^2;
      case 63, d = -sqrt(7)*cos(beta/2).^6.*sin(beta/2);
      case 64, d = cos(beta/2).^7;
    end
  case 4
    switch linidx
      case  1, d = (1+cos(beta)).^4/16;
      case  2, d = sqrt(2)/8*sin(beta).*(1+cos(beta)).^3;
      case  3, d = sqrt(7)/8*sin(beta).^2.*(1+cos(beta)).^2;
      case  4, d = sqrt(14)/8*sin(beta).^3.*(1+cos(beta));
      case  5, d = sqrt(70)/16*sin(beta).^4;
      case  6, d = sqrt(14)/8*sin(beta).^3.*(1-cos(beta));
      case  7, d = sqrt(7)/8*sin(beta).^2.*(1-cos(beta)).^2;
      case  8, d = sqrt(2)/8*sin(beta).*(1-cos(beta)).^3;
      case  9, d = (1-cos(beta)).^4/16;
      case 10, d = -sqrt(2)/8*sin(beta).*(1+cos(beta)).^3;
      case 11, d = -1/8*(1+cos(beta)).^3.*(3-4*cos(beta));
      case 12, d = -sqrt(14)/8*sin(beta).*(1+cos(beta)).^2.*(1-2*cos(beta));
      case 13, d = -sqrt(7)/8*sin(beta).^2.*(1+cos(beta)).*(1-4*cos(beta));
      case 14, d = sqrt(35)/4*sin(beta).^3.*cos(beta);
      case 15, d = sqrt(7)/8*sin(beta).^2.*(1-cos(beta)).*(1+4*cos(beta));
      case 16, d = sqrt(14)/8*sin(beta).*(1-cos(beta)).^2.*(1+2*cos(beta));
      case 17, d = 1/8*(1-cos(beta)).^3.*(3+4*cos(beta));
      case 18, d = sqrt(2)/8*sin(beta).*(1-cos(beta)).^3;
      case 19, d = sqrt(7)/8*sin(beta).^2.*(1+cos(beta)).^2;
      case 20, d = sqrt(14)/8*sin(beta).*(1+cos(beta)).^2.*(1-2*cos(beta));
      case 21, d = 1/4*(1+cos(beta)).^2.*(1-7*cos(beta)+7*cos(beta).^2);
      case 22, d = -sqrt(2)/8*sin(beta).*(1+cos(beta)).*(1+7*cos(beta)-14*cos(beta).^2);
      case 23, d = -sqrt(10)/8*sin(beta).^2.*(1-7*cos(beta).^2);
      case 24, d = -sqrt(2)/8*sin(beta).*(1-cos(beta)).*(1-7*cos(beta)-14*cos(beta).^2);
      case 25, d = 1/4*(1-cos(beta)).^2.*(1+7*cos(beta)+7*cos(beta).^2);
      case 26, d = sqrt(14)/8*sin(beta).*(1-cos(beta)).^2.*(1+2*cos(beta));
      case 27, d = sqrt(7)/8*sin(beta).^2.*(1-cos(beta)).^2;
      case 28, d = -sqrt(14)/8*sin(beta).^3.*(1+cos(beta));
      case 29, d = -sqrt(7)/8*sin(beta).^2.*(1+cos(beta)).*(1-4*cos(beta));
      case 30, d = sqrt(2)/8*sin(beta).*(1+cos(beta)).*(1+7*cos(beta)-14*cos(beta).^2);
      case 31, d = 1/8*(1+cos(beta)).*(3-6*cos(beta)-21*cos(beta).^2+28*cos(beta).^3);
      case 32, d = -sqrt(5)/4*sin(beta).*cos(beta).*(3-7*cos(beta).^2);
      case 33, d = -1/8*(1-cos(beta)).*(3+6*cos(beta)-21*cos(beta).^2-28*cos(beta).^3);
      case 34, d = -sqrt(2)/8*sin(beta).*(1-cos(beta)).*(1-7*cos(beta)-14*cos(beta).^2);
      case 35, d = sqrt(7)/8*sin(beta).^2.*(1-cos(beta)).*(1+4*cos(beta));
      case 36, d = sqrt(14)/8*sin(beta).^3.*(1-cos(beta));
      case 37, d = sqrt(70)/16*sin(beta).^4;
      case 38, d = -sqrt(35)/4*sin(beta).^3.*cos(beta);
      case 39, d = -sqrt(10)/8*sin(beta).^2.*(1-7*cos(beta).^2);
      case 40, d = sqrt(5)/4*sin(beta).*cos(beta).*(3-7*cos(beta).^2);
      case 41, d = 1/8*(3-30*cos(beta).^2+35*cos(beta).^4);
      case 42, d = -sqrt(5)/4*sin(beta).*cos(beta).*(3-7*cos(beta).^2);
      case 43, d = -sqrt(10)/8*sin(beta).^2.*(1-7*cos(beta).^2);
      case 44, d = sqrt(35)/4*sin(beta).^3.*cos(beta);
      case 45, d = sqrt(70)/16*sin(beta).^4;
      case 46, d = -sqrt(14)/8*sin(beta).^3.*(1-cos(beta));
      case 47, d = sqrt(7)/8*sin(beta).^2.*(1-cos(beta)).*(1+4*cos(beta));
      case 48, d = sqrt(2)/8*sin(beta).*(1-cos(beta)).*(1-7*cos(beta)-14*cos(beta).^2);
      case 49, d = -1/8*(1-cos(beta)).*(3+6*cos(beta)-21*cos(beta).^2-28*cos(beta).^3);
      case 50, d = sqrt(5)/4*sin(beta).*cos(beta).*(3-7*cos(beta).^2);
      case 51, d = 1/8*(1+cos(beta)).*(3-6*cos(beta)-21*cos(beta).^2+28*cos(beta).^3);
      case 52, d = -sqrt(2)/8*sin(beta).*(1+cos(beta)).*(1+7*cos(beta)-14*cos(beta).^2);
      case 53, d = -sqrt(7)/8*sin(beta).^2.*(1+cos(beta)).*(1-4*cos(beta));
      case 54, d = sqrt(14)/8*sin(beta).^3.*(1+cos(beta));
      case 55, d = sqrt(7)/8*sin(beta).^2.*(1-cos(beta)).^2;
      case 56, d = -sqrt(14)/8*sin(beta).*(1-cos(beta)).^2.*(1+2*cos(beta));
      case 57, d = 1/4*(1-cos(beta)).^2.*(1+7*cos(beta)+7*cos(beta).^2);
      case 58, d = sqrt(2)/8*sin(beta).*(1-cos(beta)).*(1-7*cos(beta)-14*cos(beta).^2);
      case 59, d = -sqrt(10)/8*sin(beta).^2.*(1-7*cos(beta).^2);
      case 60, d = sqrt(2)/8*sin(beta).*(1+cos(beta)).*(1+7*cos(beta)-14*cos(beta).^2);
      case 61, d = 1/4*(1+cos(beta)).^2.*(1-7*cos(beta)+7*cos(beta).^2);
      case 62, d = -sqrt(14)/8*sin(beta).*(1+cos(beta)).^2.*(1-2*cos(beta));
      case 63, d = sqrt(7)/8*sin(beta).^2.*(1+cos(beta)).^2;
      case 64, d = -sqrt(2)/8*sin(beta).*(1-cos(beta)).^3;
      case 65, d = 1/8*(1-cos(beta)).^3.*(3+4*cos(beta));
      case 66, d = -sqrt(14)/8*sin(beta).*(1-cos(beta)).^2.*(1+2*cos(beta));
      case 67, d = sqrt(7)/8*sin(beta).^2.*(1-cos(beta)).*(1+4*cos(beta));
      case 68, d = -sqrt(35)/4*sin(beta).^3.*cos(beta);
      case 69, d = -sqrt(7)/8*sin(beta).^2.*(1+cos(beta)).*(1-4*cos(beta));
      case 70, d = sqrt(14)/8*sin(beta).*(1+cos(beta)).^2.*(1-2*cos(beta));
      case 71, d = -1/8*(1+cos(beta)).^3.*(3-4*cos(beta));
      case 72, d = sqrt(2)/8*sin(beta).*(1+cos(beta)).^3;
      case 73, d = (1-cos(beta)).^4/16;
      case 74, d = -sqrt(2)/8*sin(beta).*(1-cos(beta)).^3;
      case 75, d = sqrt(7)/8*sin(beta).^2.*(1-cos(beta)).^2;
      case 76, d = -sqrt(14)/8*sin(beta).^3.*(1-cos(beta));
      case 77, d = sqrt(70)/16*sin(beta).^4;
      case 78, d = -sqrt(14)/8*sin(beta).^3.*(1+cos(beta));
      case 79, d = sqrt(7)/8*sin(beta).^2.*(1+cos(beta)).^2;
      case 80, d = -sqrt(2)/8*sin(beta).*(1+cos(beta)).^3;
      case 81, d = (1+cos(beta)).^4/16;
    end
end
