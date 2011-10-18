% wigner3j   Wigner 3-j symbol 
%
%   v = wigner3j(j1,j2,j3,m1,m2,m3)
%   v = wigner3j(jm1,jm2,jm3)
%   v = wigner3j(j,m)
%   v = wigner3j(jm)
%
%   Computes the value of the Wigner 3-j symbol
%
%      / j1  j2  j3 \
%      |            |
%      \ m1  m2  m3 /
%
%   Definitions for alternative input forms
%   a)  jm1 = [j1 m2], jm2 = [j2 m2], jm3 = [j3 m3]
%   b)  j = [j1 j2 j3], m = [m1 m2 m3]
%   c)  jm = [j1 j2 j3 m1 m2 m3]

function value = wigner3j(varargin)

Display = 0;

if (nargin==0), help(mfilename); return; end

% Parse input
%---------------------------------------------------
switch nargin
  case 6
    [j1,j2,j3,m1,m2,m3] = deal(varargin{:});
  case 3
    [jm1,jm2,jm3] = deal(varargin{:});
    j1 = jm1(1); m1 = jm1(2);
    j2 = jm2(1); m2 = jm2(2);
    j3 = jm3(1); m3 = jm3(2);
  case 2
    [j,m] = deal(varargin{:});
    j1 = j(1); j2 = j(2); j3 = j(3);
    m1 = m(1); m2 = m(2); m3 = m(3);
  case 1
    j = varargin{1};
    j1 = j(1); j2 = j(2); j3 = j(3);
    m1 = j(4); m2 = j(5); m3 = j(6);
  otherwise
    error('Wrong number of parameters!');
end

% Reject nonphysical parameters
%--------------------------------------------------
jjjmmm = [j1 j2 j3 m1 m2 m3];
if any(~isint(2*jjjmmm))
  error('Nonphysical parameters. All parameters should be integers or half-integers.');
end

if any([j1 j2 j3]<0)
  error('Nonphysical parameters. All j should satisfy j>=0');
end

if abs(m1)>j1 | abs(m2)>j2 | abs(m3)>j3
  error('Nonphysical parameters. m should be one of -j,-j+1,...,j-1,j.');
end

if ~isint(j1-m1) | ~isint(j2-m2) | ~isint(j3-m3)
  error('Nonphysical parameters. m should be one of -j,-j+1,...,j-1,j.');
end


% Check for zero conditions
%--------------------------------------------------
% (i) The ms must add up to zero.
if (m1+m2+m3~=0)
  value = 0;
  return;
end

% (ii) The js must satisfy the triangle relations.
if ~istriangle(j1,j2,j3)
  value = 0;
  return;
end

% Nonzero value: computation
%--------------------------------------------------

% Value for all zeros
if ~any(jjjmmm)
  value = 1;
  return;
end

% Fast formula for any j == 2
if (j2==2),
  value = fastwigner(j1,j2,j3,m1,m2,m3);
  return;
end
if (j3==2)
  value = fastwigner(j2,j3,j1,m2,m3,m1);
  return;
end
if (j1==2)
  value = fastwigner(j3,j1,j2,m3,m1,m2);
  return;
end

% Values for [4 4 4; 0 0 0] and [4 4 0; 0 0 0] and [4 0 4; 0 0 0]
mzero = ~any([m1,m2,m3]);
if mzero
  if (j1==4)
    if (j2==4)&(j3==4)
      value = sqrt(18/1001);
      return;
    elseif (j2==4)&(j3==0)
      value = 1/3;
      return;
    elseif (j2==0)&(j3==4)
      value = 1/3;
      return;
    end
  else % general routine for [j1,j2,j3;0 0 0]
    % see Tuzun (1998), p.115, Eq.(12) (typo in Eq.(13))
    % the formula from Edmonds p.125 is more prone to overflow errors
    J = j1+j2+j3;
    if mod(J,2)==1, value = 0; return; end
    CBA = sort([-j1+j2+j3,j1-j2+j3,j1+j2-j3]);
    C = CBA(1); B = CBA(2); A = CBA(3);
    v = 1/(J+1);
    for i=1:B/2, v = v*(B/2+i)*(A/2+i)^2/(A+i)/(A+B/2+i)/i; end
    for i=1:C/2, v = v*(C/2+i)*(A/2+B/2+i)^2/(A+B+i)/(A+B+C/2+i)/i; end
    value = (-1)^(J/2)*sqrt(v);
    return;
  end
end

% General computation
%==================================================================

% Formula from Eq. (1)
% Lai and Chiu, Computer Physics Communications 61 (1990) 350-360

% Pre-factor
%-----------------------------------------
PrefactorMethod = '1';
switch PrefactorMethod
  case '0' % direct evaluation of expression
    v = factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)*...
      factorial(j3+m3)*factorial(j3-m3)/factorial(j1+j2-j3)/factorial(j1-j2+j3)/...
      factorial(-j1+j2+j3)/factorial(j1+j2+j3+1);
    prefactor = v;
  case '1' % logarithmic form
    v = facln(j1+m1) + facln(j1-m1) + facln(j2+m2) + facln(j2-m2) + ...
      facln(j3+m3) + facln(j3-m3) - facln(j1+j2+j3+1) - ...
      facln(j1+j2-j3) - facln(j1-j2+j3) - facln(-j1+j2+j3);
    prefactor = exp(v/2);
  case '2' % using prime number basis
    v = mult_pr(1,fac_pr(j1+m1));
    v = mult_pr(v,fac_pr(j1-m1));
    v = mult_pr(v,fac_pr(j2+m2));
    v = mult_pr(v,fac_pr(j2-m2));
    v = mult_pr(v,fac_pr(j3+m3));
    v = mult_pr(v,fac_pr(j3-m3));
    v = divi_pr(v,fac_pr(+j1+j2-j3));
    v = divi_pr(v,fac_pr(+j1-j2+j3));
    v = divi_pr(v,fac_pr(-j1+j2+j3));
    v = divi_pr(v,fac_pr(j1+j2+j3+1));
    prefactor = pr2num(v/2);
  otherwise
    error('Unknown prefactor computation method.');
end

% Binomial sum
%-----------------------------------------
tmin = max([0,j1-j3+m2,j2-j3-m1]);
tmax = min([j1+j2-j3,j1-m1,j2+m2]);

BinsumMethod = 1;
if max([j1 j2 j3])>20
  BinsumMethod = 2;
end
switch BinsumMethod
  case 1 % logarithmic form
    binsum = 0;
    for t = tmin:tmax
      p = binoln(j1+j2-j3,t) + binoln(j1-j2+j3,j1-m1-t) + binoln(-j1+j2+j3,j2+m2-t);
      p = (-1)^t*exp(p);
      binsum = binsum + p;
    end
  case 2 % direct evaluation using Java BigInteger class
    t = tmin;
    p = binom_bi(j1+j2-j3,t);
    p = p.multiply(binom_bi( j1-j2+j3,j1-m1-t));
    p = p.multiply(binom_bi(-j1+j2+j3,j2+m2-t));
    p = mult(p,(-1)^t);
    binsum = p;
    for t = tmin+1:tmax
      p = mult(p,j1+j2-j3-t+1); p = divi(p,t);
      p = mult(p,j1-m1-t+1);    p = divi(p,-j2+j3+m1+t);
      p = mult(p,j2+m2-t+1);    p = divi(p,-j1+j3-m2+t);
      p = mult(p,-1);
      binsum = binsum.add(p);
    end
    binsum = binsum.doubleValue;
  otherwise
    error('Unknown computation method for binomial sum.');
end

value = (-1)^(j1-j2-m3)*prefactor*binsum;

%---------------
%   % Formula (29) from Tuzun 1998 (Comput.Phys.Commun.)
% 
%   alpha = [j2-j3-m1,0,j1-j3+m2];
%   beta = [j2+m2,j1+j2-j3,j1-m1];
%   A = sort(alpha); A1 = A(1); A2 = A(2); A3 = A(3);
%   B = sort(beta);  B1 = B(1); B2 = B(2); B3 = B(3);
%   
%   %  buggy...
%   
%   %p = 1;
%   %for k=1:A3-A2, p = p * (B3-A3+k)/(B2+B3-A1-A3+k); end
%   %for k=1:B2-B1, p = p * (B1-A1+k)/(B1-A2+k)*(B1-A3+k)/(B1+B3-A1-A3+k); end
%   %for k=1:B1-A3, p = p /(B3-A1+k)/(B2+B3-A1-A2+k); end
%   %p = fac(B1-A3)*sqrt(p);
%   
%   J = j1+j2+j3;
%   p = fac(A1-B1)*fac(A1-B2)/fac(A1-B3)*...
%       fac(A2-B1)*fac(A2-B2)/fac(A2-B3)*...
%       fac(A3-B1)*fac(A3-B2)/fac(A3-B3)/fac(J+1);
%   p = sqrt(p);
%   
%   t = (-1)^A3*bino(B1-A1,B1-A3)*bino(B2-A2,B2-A3);
%   s = t;
%   for k = A3+1:B1
%     t = -t*(B1+1-k)/(k-A1)*(B2+1-k)/(k-A2)*(B3+1-k)/(k-A3);
%     s = s + t;
%   end
%   s
%   
%   value = (-1)^(j1-j2-m3)*p*s;
%---------------

if (Display)
  n = numel(t);
  fprintf('sum over %d terms: %d facln in total\n',n,n*6+10);
end

% Remove numerical noise on zero results
if abs(value)<1e-10
  value = 0;
end

return

%============================================================


%===========================================================

function v = mult(p,a)
v = p.multiply(bi(a));

function v = divi(p,a)
v = p.divide(bi(a));

function v = binom_bi(n,k)
v = bi(1);
for q = n-k+1:n
  v = v.multiply(bi(q));
end
for q = 1:k
  v = v.divide(bi(q));
end

function v = bi(d)
v = java.math.BigInteger(sprintf('%d',d));

%--------------------------------------------------------------

function c = mult_pr(a,b)
if numel(a)==1, a = pr(a); end
if numel(b)==1, b = pr(b); end
c = a + b;

function c = divi_pr(a,b)
if numel(a)==1, a = pr(a); end
if numel(b)==1, b = pr(b); end
c = a - b;

function n_pr = pr(n)
v = factor(n);
persistent p; if isempty(p), p = primes(20000); end
n_pr = zeros(size(p));
if n==1, return; end
k = 1;
q = 1;
while q<=numel(v)
  while p(k)<v(q)
    k = k + 1;
  end
  while p(k)==v(q)
    n_pr(k) = n_pr(k)+1;
    q = q+1;
    if q>numel(v); break, end;
  end
end

function n = pr2num(n_pr)
persistent p; if isempty(p), p = primes(20000); end
n = prod(p.^n_pr);

function v = fac_pr(n)
v = pr(1);
for k = 2:n
  v = mult_pr(v,k);
end

%-----------------------------------------------------------------

function v = facln(x) % Logarithm of factorial
v = gammaln(x+1);

function v = binoln(n,k) % logarihtm of binomial coefficient
v = facln(n) - facln(k) - facln(n-k);

%-----------------------------------------------------------------


function w3j = fastwigner(j1,j2,j3,m1,m2,m3)
% Expressions taken from
% A.R.Edmonds, Angular Momentum, Princeton University Press, 1957
% Table 2, p.125-127

%----------------------------------------------------------------------
% Permute variables if necessary to get m2 => 0 and j1 <= j3;
% Keep track of phase
%----------------------------------------------------------------------

if (mod(j1+j2+j3,2)==0)
  parity = 1;
else
  parity = -1;
end

if (m2<0)
  m1 = -m1;
  m2 = -m2;
  m3 = -m3;
  phase = parity;
else
  phase = 1;
end

if (j1>j3)
  k = j1; j1 = j3; j3 = k;
  k = m1; m1 = m3; m3 = k;
  phase = phase*parity;
end

if (mod(j1-m3,2)~=0)
  phase = -phase;
end

%----------------------------------------------------------------------
%     Calculate Wigner 3-j symbols
%----------------------------------------------------------------------

jdelta = j3-j1;
x = 2*j1 - 1 + jdelta;
y = j1 - m3;
z = j1 + m3;

if (j2==0)

  w3j = 1/sqrt(2*j1+1);

elseif (j2==2)

  tmp2 = x*(x+1)*(x+2)*(x+3)*(x+4);
  if (m2==0)
    if (jdelta==0)
      tmp1 = 2*(3*m3*m3-j1*(j1+1));
      w3j = tmp1/sqrt(tmp2);
    elseif (jdelta==1)
      tmp1  = 6*(z+1)*(y+1);
      w3j = -2*m3*sqrt(tmp1/tmp2);
    else
      tmp1 = 6*(z+2)*(z+1)*(y+2)*(y+1);
      w3j = sqrt(tmp1/tmp2);
    end
  elseif (m2==1)
    if (jdelta==0)
      tmp1 = 6*(z+1)*y;
      w3j = (2*m3+1)*sqrt(tmp1/tmp2);
    elseif (jdelta==1)
      tmp1 = y*(y+1);
      w3j = -(2*j1+4*m3+4)*sqrt(tmp1/tmp2);
    else
      tmp1 = (z+2)*y*(y+1)*(y+2);
      w3j = 2*sqrt(tmp1/tmp2);
    end
  else
    if (jdelta==0)
      tmp1 = 6*(y-1)*y*(z+1)*(z+2);
      w3j = sqrt(tmp1/tmp2);
    elseif (jdelta==1)
      tmp1 = (y-1)*y*(y+1)*(z+2);
      w3j = -2*sqrt(tmp1/tmp2);
    else
      tmp1 = (y-1)*y*(y+1)*(y+2);
      w3j = sqrt(tmp1/tmp2);
    end
  end

else

  if (m2==0) & (jdelta==0)
    tmp2 = j1*(j1+1)*(x+2);
  else
    tmp2 = (x+1)*(x+2)*(x+3);
  end
  if (m2==0)
    if (jdelta==0)
      w3j = m3/sqrt(tmp2);
    else
      tmp1 = 2*(y+1)*(z+1);
      w3j = -sqrt(tmp1/tmp2);
    end
  else
    if (jdelta==0)
      tmp1 = 2*y*(z+1);
      w3j = sqrt(tmp1/tmp2);
    else
      tmp1 = y*(y+1);
      w3j = -sqrt(tmp1/tmp2);
    end
  end

end

w3j = w3j*phase;

return

%============================================================

function ok = isint(x)
% Checks whether parameter is integer
ok = x==floor(x);
return

function yes = istriangle(a,b,c)
% Triangle condition checker
yes = (a+b>=c) & (b+c>=a) & (c+a>=b);
return

%============================================================
% Test code

maxJ = 6;
for a = 0:maxJ
  for b = 0:maxJ
    for c = 0:maxJ
      v = wigner3j([a,b,c],[0,0,0]);
      if (v~=0)
        fprintf('%1d %1d %1d    %+.6f\n',a,b,c,v);
      end
    end
  end
end
