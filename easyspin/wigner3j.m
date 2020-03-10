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
%   c)  jjjmmm = [j1 j2 j3; m1 m2 m3]

function value = wigner3j(varargin)

if nargin==0, help(mfilename); return; end

% Parse input
%---------------------------------------------------
if ischar(varargin{end})
  Method = varargin{end};
  inputs = varargin(1:end-1);
else
  Method = '';
  inputs = varargin;
end

switch numel(inputs)
  case 6
    [j1,j2,j3,m1,m2,m3] = deal(inputs{:});
  case 3
    [jm1,jm2,jm3] = deal(inputs{:});
    j1 = jm1(1); m1 = jm1(2);
    j2 = jm2(1); m2 = jm2(2);
    j3 = jm3(1); m3 = jm3(2);
  case 2
    [j,m] = deal(inputs{:});
    j1 = j(1); j2 = j(2); j3 = j(3);
    m1 = m(1); m2 = m(2); m3 = m(3);
  case 1
    jm = inputs{1};
    if ~all(size(jm)==[2 3])
      error('If all J and M are supplied in an array, the array must be 2x3.');
    end
    j = jm(1,:);
    m = jm(2,:);
    j1 = j(1); j2 = j(2); j3 = j(3);
    m1 = m(1); m2 = m(2); m3 = m(3);
  otherwise
    error('Wrong number of parameters!');
end

if isempty(Method)
  Method = 'f+';
  if max([j1 j2 j3])>20
    Method = 'b+';
  end
end

isint = @(x) x==floor(x);
istriangle = @(a,b,c) (a+b>=c) && (b+c>=a) && (c+a>=b);

% Reject nonphysical parameters
%-------------------------------------------------------------------------------
jjjmmm = [j1 j2 j3 m1 m2 m3];
if any(~isint(2*jjjmmm))
  error('Nonphysical parameters. All parameters must be integers or half-integers.');
end

if j1<0
  error('Nonphysical parameter. j1 must satisfy j1>=0.');
end
if j2<0
  error('Nonphysical parameter. j2 must satisfy j2>=0.');
end
if j3<0
  error('Nonphysical parameter. j3 must satisfy j3>=0.');
end

if ~isint(j1-m1)
  error('Nonphysical parameter. m1 must be one of -j1,-j1+1,...,j1-1,j1.');
end
if ~isint(j2-m2)
  error('Nonphysical parameter. m2 must be one of -j2,-j2+1,...,j2-1,j2.');
end
if ~isint(j3-m3)
  error('Nonphysical parameter. m3 must be one of -j3,-j3+1,...,j3-1,j3.');
end

% Check for zero conditions
%-------------------------------------------------------------------------------
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
%-------------------------------------------------------------------------------

% Value for [0 0 0; 0 0 0]
if all(jjjmmm==0)
  value = 1;
  return
end

% Use fast explicit expressions if any j<=2
useFastExpressions = any(Method=='+');
if useFastExpressions
  if j1<=2 || j2<=2 || j3<=2
    value = fastwigner(j1,j2,j3,m1,m2,m3);
    return
  end
end

% Values for [j1 j2 j3; 0 0 0]
if useFastExpressions
  if m1==0 && m2==0 && m3==0
    % General routine for [j1 j2 j3; 0 0 0]
    % see Tuzun, Burkhardt, Secrest
    % Accurate computation of individual and tables of 3-j and 6-j symbols
    % Computer Physics Communications 112, 112-148 (1998)
    % https://doi.org/10.1016/S0010-4655(98)00065-4
    % p.115, Eq.(12) (typo in Eq.(13))
    % The expression from Edmonds p.125 is more prone to overflow errors.
    J = j1+j2+j3;
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
%===============================================================================
% Formula from Eq. (1)
% Lai and Chiu, Computer Physics Communications 61 (1990) 350-360
% https://doi.org/10.1016/0010-4655(90)90049-7

facln = @(x)gammaln(x+1); % Logarithm of factorial
binoln = @(n,k)facln(n)-facln(k)-facln(n-k); % Logarithm of binomial coefficient

tmin = max([0,j1-j3+m2,j2-j3-m1]);
tmax = min([j1+j2-j3,j1-m1,j2+m2]);
if any(Method=='f')
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

elseif any(Method=='b')
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
  
else
  
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
function val = fastwigner(j1,j2,j3,m1,m2,m3)
% Implements explicit expressions for min([j1 j2 j3])<=2.
% Expressions taken from
%   A.R.Edmonds, Angular Momentum, Princeton University Press, 1957
%   Table 2, p.125-127

phase = 1; % to keep track of overall sign

parity = (-1)^(j1+j2+j3);

% If needed, swap j3/m3 with j1/m1 or j2/m2 to get j3<=2
if j1<=2
  k = j3; j3 = j1; j1 = k;
  k = m3; m3 = m1; m1 = k;
  phase = phase*parity;
elseif j2<=2
  k = j3; j3 = j2; j2 = k;
  k = m3; m3 = m2; m2 = k;
  phase = phase*parity;
elseif j3<=2
  % ok
else
  error('At least one J must be <=2.');
end

% If needed, invert all m to get m3>=0
if m3<0
  m1 = -m1;
  m2 = -m2;
  m3 = -m3;
  phase = phase*parity;
end

% If needed, swap j1/m1 and j2/m2 to get j1>=j2
if j1<j2
  k = j1; j1 = j2; j2 = k;
  k = m1; m1 = m2; m2 = k;
  phase = phase*parity;
end

% Calculate Wigner 3-j symbols
J = j2;
M = m1;
JmM = J - M;
JpM = J + M;
jdelta = j1-j2;
x = j1 + j2 - 1;

phase = phase*(-1)^(J-M);

if j3==0

  val = 1/sqrt(2*j1+1);
    
elseif j3==2

  tmp2 = x*(x+1)*(x+2)*(x+3)*(x+4);
  if m3==0
    if jdelta==0
      val = 2*(3*M^2-J*(J+1))/sqrt(tmp2);
    elseif jdelta==1
      tmp1  = 6*(JpM+1)*(JmM+1);
      val = -2*M*sqrt(tmp1/tmp2);
    else
      tmp1 = 6*(JpM+2)*(JpM+1)*(JmM+2)*(JmM+1);
      val = sqrt(tmp1/tmp2);
    end
  elseif m3==1
    if jdelta==0
      tmp1 = 6*(JpM+1)*JmM;
      val = (1+2*M)*sqrt(tmp1/tmp2);
    elseif jdelta==1
      tmp1 = (JmM+1)*JmM;
      val = -2*(J+2*M+2)*sqrt(tmp1/tmp2);
    else % jdelta==2
      tmp1 = (JpM+2)*(JmM+2)*(JmM+1)*JmM;
      val = 2*sqrt(tmp1/tmp2);
    end
  else % jdelta==2
    if jdelta==0
      tmp1 = 6*(JmM-1)*JmM*(JpM+1)*(JpM+2);
      val = sqrt(tmp1/tmp2);
    elseif jdelta==1
      tmp1 = (JmM-1)*JmM*(JmM+1)*(JpM+2);
      val = -2*sqrt(tmp1/tmp2);
    else % jdelta==2
      tmp1 = (JmM-1)*JmM*(JmM+1)*(JmM+2);
      val = sqrt(tmp1/tmp2);
    end
  end

elseif j3==1/2

  val = -1i*sqrt((JmM+1/2)/(2*J+2)/(2*J+1));

elseif j3==1

  if m3==0 && jdelta==0
    tmp2 = J*(J+1)*(x+2);
  else
    tmp2 = (x+1)*(x+2)*(x+3);
  end
  if m3==0
    if jdelta==0
      val = M/sqrt(tmp2);
    else % jdelta==1
      tmp1 = (JmM+1)*(JpM+1)*2;
      val = -sqrt(tmp1/tmp2);
    end
  else % m3==1
    if jdelta==0
      tmp1 = JmM*(JpM+1)*2;
      val = sqrt(tmp1/tmp2);
    else % jdelta==1
      tmp1 = JmM*(JmM+1);
      val = -sqrt(tmp1/tmp2);
    end
  end
  
elseif j3==3/2
  
  if m3==3/2
    if jdelta==3/2
      val = 1i*sqrt((JmM-1/2)*(JmM+1/2)*(JmM+3/2)/(2*J+4)/(2*J+3)/(2*J+2)/(2*J+1));
    else % jdelta==1/2
      val = -1i*sqrt(3*(JmM-1/2)*(JmM+1/2)*(JpM+3/2)/(2*J+3)/(2*J+2)/(2*J+1)/(2*J));
    end
  else % m3==1/2
    if jdelta==3/2
      val = 1i*sqrt(3*(JmM+1/2)*(JmM+3/2)*(JpM+3/2)/(2*J+4)/(2*J+3)/(2*J+2)/(2*J+1));
    else % jdelta==1/2
      val = -1i*sqrt((JmM+1/2)/(2*J+3)/(2*J+2)/(2*J+1)/(2*J))*(J+3*M+3/2);
    end
  end
  
else
  
  error('j3 must be 0, 1/2, 1, 3/2, or 2.');
  
end

% Apply overall phase
val = val*phase;

return
