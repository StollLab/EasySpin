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

if nargin==0, help(mfilename); return; end

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

if isempty(Method)
  Method = 1;
  if max([j1 j2 j3])>20
    Method = 2;
  end
end


isint = @(x) x==floor(x);
istriangle = @(a,b,c) (a+b>=c) && (b+c>=a) && (c+a>=b);

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
if m1+m2+m3~=0
  value = 0;
  return
end

% (ii) The js must satisfy the triangle relations.
if ~istriangle(j1,j2,j3)
  value = 0;
  return
end

% (iii) The sum of j must be even if the sum of m is zero
if m1==0 && m2==0 && m3==0 && mod(j1+j2+j3,2)
  value = 0;
  return
end

% Nonzero value: computation
%--------------------------------------------------

% Value for all zeros
if all(jjjmmm==0)
  value = 1;
  return
end

% Use fast explicit expressions for any j == 2
useFastExpressions = Method>0;
if useFastExpressions
  if j2==2
    value = fastwigner(j1,j2,j3,m1,m2,m3);
    return
  end
  if j3==2
    value = fastwigner(j2,j3,j1,m2,m3,m1);
    return
  end
  if j1==2
    value = fastwigner(j3,j1,j2,m3,m1,m2);
    return
  end
end

% Values for [j1,j2,j3;0 0 0]
if m1==0 && m2==0 && m3==0
  J = j1+j2+j3;
  if mod(J,2)
    value = 0;
    return
  end
  % Values for [4 4 4; 0 0 0] and [4 4 0; 0 0 0] and [4 0 4; 0 0 0]
  if j1==4
    if j2==4 && j3==4
      value = sqrt(18/1001);
      return
    elseif j2==4 && j3==0
      value = 1/3;
      return
    elseif j2==0 && j3==4
      value = 1/3;
      return
    end
  else
    % General routine for [j1,j2,j3;0 0 0]
    % see Tuzun, Burkhardt, Secrest
    % Accurate computation of individual and tables of 3-j and 6-j symbols
    % Computer Physics Communications 112, 112-148 (1998)
    % https://doi.org/10.1016/S0010-4655(98)00065-4
    % p.115, Eq.(12) (typo in Eq.(13))
    % The formula from Edmonds p.125 is more prone to overflow errors.
    if J<100000
      CBA = sort([-j1+j2+j3,j1-j2+j3,j1+j2-j3]);
      C = CBA(1); B = CBA(2); A = CBA(3);
      v = 1/(J+1);
      for i = 1:B/2
        v = v/i*(B/2+i)*(A/2+i)^2/(A+B/2+i)/(A+i);
      end
      for i = 1:C/2
        v = v/i*(C/2+i)*(A/2+B/2+i)^2/(A+B+C/2+i)/(A+B+i);
      end
      value = (-1)^(J/2)*sqrt(v);
      return
    end
  end
end

% General computation
%==================================================================
% Formula from Eq. (1)
% Lai and Chiu, Computer Physics Communications 61 (1990) 350-360
% https://doi.org/10.1016/0010-4655(90)90049-7

facln = @(x)gammaln(x+1); % Logarithm of factorial
binoln = @(n,k)facln(n)-facln(k)-facln(n-k); % Logarithm of binomial coefficient

tmin = max([0,j1-j3+m2,j2-j3-m1]);
tmax = min([j1+j2-j3,j1-m1,j2+m2]);
switch abs(Method)
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
    % sum: binomials, using Java class BigInteger/BigDecimal
    
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
    % don't merge the following three lines - Matlab 7.5 throws an error
    b = java.math.BigDecimal(binsum).movePointLeft(n).doubleValue;
    %b = b.movePointLeft(n);
    %b = b.doubleValue;
    
    prefactor_ln = ...
      facln(j1+m1) + facln(j1-m1) + facln(j2+m2) + facln(j2-m2) + ...
      facln(j3+m3) + facln(j3-m3) - facln(j1+j2+j3+1) - ...
      facln(j1+j2-j3) - facln(j1-j2+j3) - facln(-j1+j2+j3);
    
    value = (-1)^(j1-j2-m3)*exp(prefactor_ln/2+n*log(10))*b;
    
  otherwise
    error('Unknown computation method.');
end

return

%===============================================================================
%===============================================================================


%-------------------------------------------------------------------------------
function v = bi(d)
v = java.math.BigInteger(sprintf('%d',d));

%-------------------------------------------------------------------------------
function v = binom_bi(n,k)
v = bi(1);
for q = n-k+1:n
  v = v.multiply(bi(q));
end
for q = 1:k
  v = v.divide(bi(q));
end

%-------------------------------------------------------------------------------
function w3j = fastwigner(j1,j2,j3,m1,m2,m3)
% Expressions taken from
% A.R.Edmonds, Angular Momentum, Princeton University Press, 1957
% Table 2, p.125-127

%----------------------------------------------------------------------
% Permute variables if necessary to get m2 => 0 and j1 <= j3;
% Keep track of phase
%----------------------------------------------------------------------

if mod(j1+j2+j3,2)==0
  parity = 1;
else
  parity = -1;
end

if m2<0
  m1 = -m1;
  m2 = -m2;
  m3 = -m3;
  phase = parity;
else
  phase = 1;
end

if j1>j3
  k = j1; j1 = j3; j3 = k;
  k = m1; m1 = m3; m3 = k;
  phase = phase*parity;
end

if mod(j1-m3,2)~=0
  phase = -phase;
end

%----------------------------------------------------------------------
% Calculate Wigner 3-j symbols
%----------------------------------------------------------------------

jdelta = j3-j1;
x = 2*j1 - 1 + jdelta;
y = j1 - m3;
z = j1 + m3;

if j2==0

  w3j = 1/sqrt(2*j1+1);

elseif j2==2

  tmp2 = x*(x+1)*(x+2)*(x+3)*(x+4);
  if m2==0
    if jdelta==0
      tmp1 = 2*(3*m3*m3-j1*(j1+1));
      w3j = tmp1/sqrt(tmp2);
    elseif jdelta==1
      tmp1  = 6*(z+1)*(y+1);
      w3j = -2*m3*sqrt(tmp1/tmp2);
    else
      tmp1 = 6*(z+2)*(z+1)*(y+2)*(y+1);
      w3j = sqrt(tmp1/tmp2);
    end
  elseif m2==1
    if jdelta==0
      tmp1 = 6*(z+1)*y;
      w3j = (2*m3+1)*sqrt(tmp1/tmp2);
    elseif jdelta==1
      tmp1 = y*(y+1);
      w3j = -(2*j1+4*m3+4)*sqrt(tmp1/tmp2);
    else
      tmp1 = (z+2)*y*(y+1)*(y+2);
      w3j = 2*sqrt(tmp1/tmp2);
    end
  else
    if jdelta==0
      tmp1 = 6*(y-1)*y*(z+1)*(z+2);
      w3j = sqrt(tmp1/tmp2);
    elseif jdelta==1
      tmp1 = (y-1)*y*(y+1)*(z+2);
      w3j = -2*sqrt(tmp1/tmp2);
    else
      tmp1 = (y-1)*y*(y+1)*(y+2);
      w3j = sqrt(tmp1/tmp2);
    end
  end

else

  if m2==0 && jdelta==0
    tmp2 = j1*(j1+1)*(x+2);
  else
    tmp2 = (x+1)*(x+2)*(x+3);
  end
  if m2==0
    if jdelta==0
      w3j = m3/sqrt(tmp2);
    else
      tmp1 = 2*(y+1)*(z+1);
      w3j = -sqrt(tmp1/tmp2);
    end
  else
    if jdelta==0
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
