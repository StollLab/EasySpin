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
%   a)  jm1 = [j1 m1], jm2 = [j2 m2], jm3 = [j3 m3]
%   b)  jjj = [j1 j2 j3], mmm = [m1 m2 m3]
%   c)  jjjmmm = [j1 j2 j3; m1 m2 m3]

function value = wigner3j(varargin)

if nargin==0, help(mfilename); return; end

% Parse input
%---------------------------------------------------
% Undocumented: last input 'e' (default) uses explicit expressions where
% available, '' uses the recursion for all cases.
% (Written with scalar operations only, since this part dominates the run
% time for most calls.)
nInputs = nargin;
if ischar(varargin{end})
  Method = varargin{end};
  nInputs = nInputs-1;
else
  Method = 'e';
end
if ~isempty(Method) && ~strcmp(Method,'e')
  error('Unknown method ''%s''.',Method);
end
useExplicitExpressions = ~isempty(Method);

switch nInputs
  case 6
    j1 = varargin{1}; j2 = varargin{2}; j3 = varargin{3};
    m1 = varargin{4}; m2 = varargin{5}; m3 = varargin{6};
  case 3
    jm1 = varargin{1}; jm2 = varargin{2}; jm3 = varargin{3};
    j1 = jm1(1); m1 = jm1(2);
    j2 = jm2(1); m2 = jm2(2);
    j3 = jm3(1); m3 = jm3(2);
  case 2
    j = varargin{1}; m = varargin{2};
    j1 = j(1); j2 = j(2); j3 = j(3);
    m1 = m(1); m2 = m(2); m3 = m(3);
  case 1
    jm = varargin{1};
    if ~isequal(size(jm),[2 3])
      error('If all J and M are supplied in an array, the array must be 2x3.');
    end
    j1 = jm(1,1); j2 = jm(1,2); j3 = jm(1,3);
    m1 = jm(2,1); m2 = jm(2,2); m3 = jm(2,3);
  otherwise
    error('Wrong number of parameters!');
end

% Reject nonphysical parameters
%-------------------------------------------------------------------------------
if 2*j1~=floor(2*j1) || 2*j2~=floor(2*j2) || 2*j3~=floor(2*j3) || ...
   2*m1~=floor(2*m1) || 2*m2~=floor(2*m2) || 2*m3~=floor(2*m3)
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

if j1-m1~=floor(j1-m1)
  error('Nonphysical parameter. m1 must be one of -j1,-j1+1,...,j1-1,j1.');
end
if j2-m2~=floor(j2-m2)
  error('Nonphysical parameter. m2 must be one of -j2,-j2+1,...,j2-1,j2.');
end
if j3-m3~=floor(j3-m3)
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
if j1+j2<j3 || j2+j3<j1 || j3+j1<j2
  value = 0;
  return
end

% (iii) Each m must satisfy |m|<=j.
if abs(m1)>j1 || abs(m2)>j2 || abs(m3)>j3
  value = 0;
  return
end

% (iv) The sum of j must be even if all m are zero
if m1==0 && m2==0 && m3==0 && mod(j1+j2+j3,2)
  value = 0;
  return
end

% Nonzero value: computation
%-------------------------------------------------------------------------------

% Value for [0 0 0; 0 0 0]
if j1==0 && j2==0 && j3==0
  value = 1;
  return
end

% Use explicit expressions if any j<=2
if useExplicitExpressions
  if j1<=2 || j2<=2 || j3<=2
    value = wigner3j_explicit(j1,j2,j3,m1,m2,m3);
    return
  end
end

% Values for [j1 j2 j3; 0 0 0]
if useExplicitExpressions
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

% General computation: three-term recursion in j1
%===============================================================================
value = wigner3j_recursion(j1,j2,j3,m1,m2,m3);

end

%===============================================================================
%===============================================================================


%-------------------------------------------------------------------------------
function value = wigner3j_recursion(j1,j2,j3,m1,m2,m3)
% Computes the 3j symbol via the three-term recursion in j1 of
%   K. Schulten, R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
%   https://doi.org/10.1063/1.522426
% All 3j symbols with j1 = j1min..j1max are computed, with forward recursion
% from j1min and backward recursion from j1max, joined in the classically
% allowed region. Normalization and sign are fixed via
%   sum_j1 (2*j1+1)*(3j symbol)^2 = 1
%   sign of 3j symbol for j1max = (-1)^(j2-j3-m1)
% Written with scalar loops only, since vector operations on the short
% sequences involved are slower due to their overhead.

% Recurse over the largest j, since this gives the shortest recursion.
% An odd permutation of columns gives a phase factor (-1)^(j1+j2+j3).
phase = 1;
if j2>j1 && j2>=j3
  tmp = j1; j1 = j2; j2 = tmp;
  tmp = m1; m1 = m2; m2 = tmp;
  phase = (-1)^(j1+j2+j3);
elseif j3>j1
  tmp = j1; j1 = j3; j3 = tmp;
  tmp = m1; m1 = m3; m3 = tmp;
  phase = (-1)^(j1+j2+j3);
end

jmin = max(abs(j2-j3),abs(m1));
jmax = j2+j3;
N = jmax-jmin+1;
iTarget = j1-jmin+1;

if N==1
  value = phase*(-1)^(j2-j3-m1)/sqrt(2*j1+1);
  return
end

% Recursion: cUp(j)*f(j+1) + cMid(j)*f(j) + cDown(j)*f(j-1) = 0, with
%   cUp(j) = j*A(j+1), cDown(j) = (j+1)*A(j)
%   A(j) = sqrt((j^2-(j2-j3)^2)*((j2+j3+1)^2-j^2)*(j^2-m1^2))
%   cMid(j) = -(2j+1)*(c1 - j(j+1)*dm)
d2 = (j2-j3)^2;
s2 = (j2+j3+1)^2;
mm = m1^2;
c1 = (j2*(j2+1)-j3*(j3+1))*m1;
dm = m3-m2;

bigValue = 1e100; % rescaling threshold to avoid overflow

% Forward recursion from jmin, until the values stop increasing (i.e.
% until past the classically forbidden region at small j1)
f = zeros(1,N);
f(1) = 1;
j = jmin;
jj = (j+1)^2;
Ajp = sqrt(max(0,(jj-d2)*(s2-jj)*(jj-mm))); % A(jmin+1)
if jmin==0
  % j2==j3 and m1==0: first recursion equation is trivial, use explicit
  % values for j1 = 0 and 1 instead
  f(2) = m2/sqrt(j2*(j2+1));
else
  f(2) = (2*j+1)*(c1-j*(j+1)*dm)/(j*Ajp);
end
k = N;
for i = 2:N-1
  j = jmin+i-1;
  Aj = Ajp;
  jj = (j+1)^2;
  Ajp = sqrt(max(0,(jj-d2)*(s2-jj)*(jj-mm)));
  f(i+1) = ((2*j+1)*(c1-j*(j+1)*dm)*f(i) - (j+1)*Aj*f(i-1))/(j*Ajp);
  if abs(f(i+1))>bigValue
    f(1:i+1) = f(1:i+1)/bigValue;
  end
  if abs(f(i+1))<abs(f(i))
    k = i;
    break
  end
end

if k==N
  % Forward recursion reached jmax: normalize and fix sign
  nrm = 0;
  for i = 1:N
    nrm = nrm + (2*(jmin+i-1)+1)*f(i)^2;
  end
  value = f(iTarget)/sqrt(nrm);
  if (f(N)>0) ~= ((-1)^(j2-j3-m1)>0)
    value = -value;
  end
  value = phase*value;
  return
end

% Backward recursion from jmax down to k-1
b = zeros(1,N);
b(N) = 1;
j = jmax;
jj = j^2;
Aj = sqrt(max(0,(jj-d2)*(s2-jj)*(jj-mm))); % A(jmax)
b(N-1) = (2*j+1)*(c1-j*(j+1)*dm)/((j+1)*Aj);
for i = N-1:-1:max(k,2)
  j = jmin+i-1;
  Ajp = Aj;
  jj = j^2;
  Aj = sqrt(max(0,(jj-d2)*(s2-jj)*(jj-mm)));
  b(i-1) = ((2*j+1)*(c1-j*(j+1)*dm)*b(i) - j*Ajp*b(i+1))/((j+1)*Aj);
  if abs(b(i-1))>bigValue
    b(i-1:N) = b(i-1:N)/bigValue;
  end
end

% Join: scale backward values onto forward values via least-squares fit over
% the overlap k-1..k+1
num = 0;
den = 0;
for i = max(k-1,1):k+1
  num = num + f(i)*b(i);
  den = den + b(i)^2;
end
scale = num/den;

% Normalize (values 1..k from forward, k+1..N from scaled backward recursion)
nrm = 0;
for i = 1:k
  nrm = nrm + (2*(jmin+i-1)+1)*f(i)^2;
end
nrmB = 0;
for i = k+1:N
  nrmB = nrmB + (2*(jmin+i-1)+1)*b(i)^2;
end
nrm = nrm + scale^2*nrmB;

if iTarget<=k
  value = f(iTarget);
else
  value = scale*b(iTarget);
end
value = value/sqrt(nrm);
% sign of the 3j symbol for jmax is (-1)^(j2-j3-m1); b(N) = 1
if (scale>0) ~= ((-1)^(j2-j3-m1)>0)
  value = -value;
end
value = phase*value;

end

%-------------------------------------------------------------------------------
function val = wigner3j_explicit(j1,j2,j3,m1,m2,m3)
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
  else % m3==2
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

end
