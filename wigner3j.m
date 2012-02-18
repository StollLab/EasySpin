% wigner3j   Wigner 3-j symbol 
%
%   v = wigner3j(j1,j2,j3,m1,m2,m3)
%   v = wigner3j(jm1,jm2,jm3)
%   v = wigner3j(jjj,mmm)
%   v = wigner3j(jjjmmm)
%
%   Computes the value of the Wigner 3-j symbol
%
%      / j1  j2  j3 \
%      |            |
%      \ m1  m2  m3 /
%
%   Definitions for alternative input forms
%   a)  jm1 = [j1 m2], jm2 = [j2 m2], jm3 = [j3 m3]
%   b)  jjj = [j1 j2 j3], mmm = [m1 m2 m3]
%   c)  jjjmmm = [j1 j2 j3 m1 m2 m3]

function value = wigner3j(varargin)

if (nargin==0), help(mfilename); return; end

% Parse input
%---------------------------------------------------
Method = [];
switch nargin
  case 7
    [j1,j2,j3,m1,m2,m3,Method] = deal(varargin{:});
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

if abs(m1)>j1 || abs(m2)>j2 || abs(m3)>j3
  error('Nonphysical parameters. m should be one of -j,-j+1,...,j-1,j.');
end

if ~isint(j1-m1) || ~isint(j2-m2) || ~isint(j3-m3)
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
if all(jjjmmm==0)
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
mzero = all([m1,m2,m3]==0);
if mzero
  if (j1==4)
    if (j2==4)&&(j3==4)
      value = sqrt(18/1001);
      return;
    elseif (j2==4)&&(j3==0)
      value = 1/3;
      return;
    elseif (j2==0)&&(j3==4)
      value = 1/3;
      return;
    end
  else % general routine for [j1,j2,j3;0 0 0]
    % see Tuzun (1998), p.115, Eq.(12) (typo in Eq.(13))
    % the formula from Edmonds p.125 is more prone to overflow errors
    J = j1+j2+j3;
    if mod(J,2)==1, value = 0; return; end
    if J<100000
      CBA = sort([-j1+j2+j3,j1-j2+j3,j1+j2-j3]);
      C = CBA(1); B = CBA(2); A = CBA(3);
      v = 1/(J+1);
      for i=1:B/2, v = v/i*(B/2+i)*(A/2+i)^2/(A+B/2+i)/(A+i); end
      for i=1:C/2, v = v/i*(C/2+i)*(A/2+B/2+i)^2/(A+B+C/2+i)/(A+B+i); end
      value = (-1)^(J/2)*sqrt(v);
      return
    end
  end
end

% General computation
%==================================================================
% Formula from Eq. (1)
% Lai and Chiu, Computer Physics Communications 61 (1990) 350-360

if isempty(Method)
  Method = 1;
  if max([j1 j2 j3])>20
    Method = 2;
  end
end

tmin = max([0,j1-j3+m2,j2-j3-m1]);
tmax = min([j1+j2-j3,j1-m1,j2+m2]);
switch Method
  case 1
    % prefactor: logarithmic
    % sum: each term logarithmic
    v = facln(j1+m1) + facln(j1-m1) + facln(j2+m2) + facln(j2-m2) + ...
      facln(j3+m3) + facln(j3-m3) - facln(j1+j2+j3+1) - ...
      facln(j1+j2-j3) - facln(j1-j2+j3) - facln(-j1+j2+j3);
    binsum = 0;
    for t = tmin:tmax
      p = binoln(j1+j2-j3,t) + binoln(j1-j2+j3,j1-m1-t) + binoln(-j1+j2+j3,j2+m2-t);
      p = (-1)^t*exp(p+v/2);
      binsum = binsum + p;
    end
    value = (-1)^(j1-j2-m3)*binsum;

  case 2
    % prefactor: logarithmic
    % sum: binomials, using Java BigInteger/BigDecimal
    for q=j1+j2+j3:-1:1
      %bigint(q) = java.math.BigInteger(sprintf('%d',q));
    end
    
    % binomial sum
    t = tmin;
    p = binom_bi(j1+j2-j3,t);
    p = p.multiply(binom_bi( j1-j2+j3,j1-m1-t));
    p = p.multiply(binom_bi(-j1+j2+j3,j2+m2-t));
    p = p.multiply(bi((-1)^t));
    binsum = p;
    for t = tmin+1:tmax
      q1 = (j1+j2-j3-t+1)*(j1-m1-t+1)*(j2+m2-t+1);
      q2 = t*(-j2+j3+m1+t)*(-j1+j3-m2+t);
      p = p.multiply(bi(q1));
      p = p.divide(bi(q2));
      p = p.multiply(bi(-1));
      binsum = binsum.add(p);
    end
    n = length(binsum.toString)-1; % 10-base exponent
    b = java.math.BigDecimal(binsum).movePointLeft(n).doubleValue;
    
    prefactor_ln = ...
      facln(j1+m1) + facln(j1-m1) + facln(j2+m2) + facln(j2-m2) + ...
      facln(j3+m3) + facln(j3-m3) - facln(j1+j2+j3+1) - ...
      facln(j1+j2-j3) - facln(j1-j2+j3) - facln(-j1+j2+j3);
    
    value = (-1)^(j1-j2-m3)*exp(prefactor_ln/2+n*log(10))*b;
    
  case 3
    % prefactor: logarithmic
    % sum: binomials, using my own arbitrary-precision integers
    % (currently buggy)
    
    % binomial sum
    t = tmin;
    p = hpi_binom(j1+j2-j3,t);
    p = hpi_multiply(p,hpi_binom( j1-j2+j3,j1-m1-t));
    p = hpi_multiply(p,hpi_binom(-j1+j2+j3,j2+m2-t));
    p = hpi_multiply(p,(-1)^t);
    binsum = p;
        
    for t = tmin+1:tmax
      q1 = (j1+j2-j3-t+1)*(j1-m1-t+1)*(j2+m2-t+1);
      q2 = t*(-j2+j3+m1+t)*(-j1+j3-m2+t);
      p = hpi_multiply(p,q1);
      p = hpi_divide(p,q2);
      p = hpi_multiply(p,-1);
      binsum = hpi_add(binsum,p);
    end    
    sgn = binsum.sign;
    binsum.sign = 1;
    [b,n] = hpi_mantissaexponent(binsum);
    
    prefactor_ln = ...
      facln(j1+m1) + facln(j1-m1) + facln(j2+m2) + facln(j2-m2) + ...
      facln(j3+m3) + facln(j3-m3) - facln(j1+j2+j3+1) - ...
      facln(j1+j2-j3) - facln(j1-j2+j3) - facln(-j1+j2+j3);
    
    value = sgn*(-1)^(j1-j2-m3)*exp(prefactor_ln/2+n*log(10))*b;

  case 4
    % prefactor: logarithmic
    % sum: binomials, using John D'Errico's vpi
    
    % binomial sum
    t = vpi(tmin);
    p = nchoosek( j1+j2-j3,t)*...
        nchoosek( j1-j2+j3,j1-m1-t)*...
        nchoosek(-j1+j2+j3,j2+m2-t);
    p = p*(-1)^t;
    binsum = p;
    for t = tmin+1:tmax
      q1 = (j1+j2-j3-t+1)*(j1-m1-t+1)*(j2+m2-t+1);
      q2 = t*(-j2+j3+m1+t)*(-j1+j3-m2+t);
      p = -p*q1/q2;
      binsum = binsum + p;
    end
    
    binsum_ln = log(abs(binsum));
    
    prefactor_ln = ...
      facln(j1+m1) + facln(j1-m1) + facln(j2+m2) + facln(j2-m2) + ...
      facln(j3+m3) + facln(j3-m3) - facln(j1+j2+j3+1) - ...
      facln(j1+j2-j3) - facln(j1-j2+j3) - facln(-j1+j2+j3);
    
    value = sign(binsum)*(-1)^(j1-j2-m3)*exp(prefactor_ln/2+binsum_ln);

  otherwise
    error('Unknown computation method.');
end

return

%============================================================


%===========================================================

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

function v = binoln(n,k) % Logarithm of binomial coefficient
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

  if (m2==0) && (jdelta==0)
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


%%==============================================================
function [mantissa,exponent] = hpi_mantissaexponent(N)
if ~isstruct(N)
  if isinf(N), error('N is Inf.'); end
  N = hpi(N);
end
if N.sign==-1, error('N must be non-negative.'); end
mantissa = 0;
exponent = 0;
nDigits = numel(N.digits);
for k = 1:nDigits
  mantissa = mantissa/1000 + N.digits(k);
  exponent = exponent + 3;
end
exponent = exponent - 3;
if mantissa>=100, mantissa=mantissa/100; exponent = exponent+2; end
if mantissa>=10, mantissa=mantissa/10; exponent = exponent+1; end

function c = hpi_add(a,b)

% convert both to hpi
if ~isstruct(a), a = hpi(a); end
if ~isstruct(b), b = hpi(b); end

% bring both to the same number of digits
na = numel(a.digits);
nb = numel(b.digits);
if (na>nb), b.digits(na) = 0; end
if (nb>na), a.digits(nb) = 0; end

% add, including the sign
result = a.sign*a.digits + b.sign*b.digits;

% determine the sign of the result
hsd = find(result,1,'last');
if isempty(hsd)
  c.sign = 1;
  c.digits = 0;
  return
else
  c.sign = sign(result(hsd));
end
if (c.sign==-1), result = -result; end

base = 1000;

% carry-over
K = find((result<0) | (result>=base));
while ~isempty(K)
  if K(end) == numel(result), result(end+1) = 0; end
  olddigits = result(K);
  newdigits = mod(olddigits,base);
  result(K) = newdigits;
  carry = (olddigits - newdigits)/base;
  K = K + 1;
  result(K) = result(K) + carry;
  K = K((result(K)<0) | (result(K)>=base));
end

% trim empty digits
hsd = find(result,1,'last');
if isempty(hsd)
  c.digits = 0;
elseif length(result)>hsd
  c.digits = result(1:hsd);
else
  c.digits = result;
end

function c = hpi_divide(a,b)
c.sign = a.sign;
c.digits = [];
if (b<0), c.sign = -c.sign; end
cv = a.digits;
nDigits = find(a.digits,1,'last');
for d=nDigits:-1:1
  result(d) = fix(cv(d)/b);
  remainder = cv(d) - b*result(d);
  if (d>1)
    cv(d-1) = cv(d-1) + remainder*1000;
  end
end
c.digits = result;
return

function value = hpi(a)
if ischar(a)
  value.sign = +1;
  if a(1)=='+', value.sign = +1; a(1) = []; end
  if a(1)=='-', value.sign = -1; a(1) = []; end
  d = rem(length(a),3);
  if (d==1), a = ['00' a]; end
  if (d==2), a = ['0' a]; end
  k = 1;
  value.digits = [];
  for idx = length(a)-2:-3:1
    value.digits(k) = str2double(a(idx:idx+2));
    k = k+1;
  end  
else
  if (a<0), value.sign = -1; else value.sign = +1; end
  a = abs(a);
  k = 1;
  value.digits(1) = 0;
  while (a>0)
    value.digits(k) = rem(a,1000);
    a = fix(a/1000);
    k = k+1;
  end
end

function c = hpi_multiply(a,b)
if ~isstruct(a), a = hpi(a); end
if ~isstruct(b), b = hpi(b); end

c.sign = a.sign*b.sign;
cv = conv(a.digits,b.digits);

% carry over
nDigits = find(cv,1,'last');
cv(nDigits+1) = 0;
for i=1:nDigits
  if cv(i)>1000
    cv(i+1) = cv(i+1)+fix(cv(i)/1000);
    cv(i) = rem(cv(i),1000);
  end
end
if cv(end)==0, cv(end) = []; end
c.digits = cv;

function v = hpi_binom(n,k)
v = hpi(1);
for q = n-k+1:n, v = hpi_multiply(v,q); end
for q = 1:k, v = hpi_divide(v,q); end

function d = hpi_double(v)
nDigits = find(v.digits,1,'last');
d = 0;
for i=nDigits:-1:1
  d = d*1000 + v.digits(i);
end
d = d*v.sign;

function hpi_show(v)
fprintf([hpi_string(v) '\n']);

function str = hpi_string(v)
str = '';
if v.sign<1, str = [str '-']; end
dig = v.digits;
hsd = find(dig,1,'last');
for d=hsd:-1:1
  str = [str sprintf('%03d',dig(d))];
end
if str(1)=='0'; str = str(2:end); end
if str(1)=='0'; str = str(2:end); end
