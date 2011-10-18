% wigner6j   Wigner 6j symbol 
%
%    v = wigner6j(j1,j2,j3,J1,J2,J3)
%    v = wigner6j(jJ1,jJ2,jJ3)
%    v = wigner6j(j,J)
%    v = wigner6j(jJ)
%
%    Compute the Wigner 6j symbol.
%
%       / j1  j2  j3 \
%       \ J1  J2  J3 /  
%
%    Definitions for alternative input forms:
%    -  jJ1 = [j1 J1] etc.
%    -  j = [j1 j2 j3], J = [J1 J2 J3]
%    -  jJ = [j1 j2 j3 J1 J2 J3]

% Stefan Stoll, 31 May 2007

% Reference:
%  R.E.Tuzun, P.Burkhardt,D.Secrest
%  Accurate compuation of individual and tables of 3-j and 6-j symbols
%  Computer Physics Communications 112 (1998) 112-148
%  and references therein

function value = wigner6j(varargin)

% Parse input
%---------------------------------------------------
if (nargin==0), help(mfilename); return; end

switch nargin
  case 6
    [j1,j2,j3,j4,j5,j6] = deal(varargin{:});
  case 3
    [jJ1,jJ2,jJ3] = deal(varargin{:});
    j1 = jJ1(1); j4 = jJ1(2);
    j2 = jJ2(1); j5 = jJ2(2);
    j3 = jJ3(1); j6 = jJ3(2);
  case 2
    [j,J] = deal(varargin{:});
    j1 = j(1); j2 = j(2); j3 = j(3);
    j4 = J(1); j5 = J(2); j6 = J(3);
  case 1
    j = varargin{1};
    if numel(j)==6
      j1 = j(1); j2 = j(2); j3 = j(3);
      j4 = j(4); j5 = j(5); j6 = j(6);
    else
      if size(j,2)==6
        for k=1:size(j,1)
          value(k) = wigner6j(j(k,:));
        end
        return;
      else
        error('Wrong size of input array, must be n x 6.');
      end
    end
  otherwise
    error('Wrong number of parameters!');
end

% Reject nonphysical parameters
%--------------------------------------------------
if any(~isint(2*[j1 j2 j3 j4 j5 j6]))
 error('All parameters must be integers or half-integers.');
end

% Check for zero conditions
%--------------------------------------------------

%[j1 j2 j3; j4 j5 j6]

% A 6j symbol must satisfy 4 triangle relations to be nonzero.
if ~(istriangle(j1,j2,j3) & istriangle(j1,j5,j6) & istriangle(j4,j2,j6) & istriangle(j4,j5,j3))
  value = 0;
  return;
end

% Nonzero value: computation
%--------------------------------------------------

alpha = [j1+j2+j4+j5, j1+j3+j4+j6, j2+j3+j5+j6];
beta = [j1+j2+j3, j1+j5+j6, j2+j4+j6, j3+j4+j5];

% All methods implement Formula (80) on p.131 of Tuzun 1998.
% Method 0: logarithmic form
% Method 1: as is
% Method 2: common factor pulled out of sum, recursive computation
%  of terms in the sum, prefactor rewritten in terms of A and B

%ComputationMethod = 2;
%if (ComputationMethod==0)
%
%  n = max(beta):min(alpha);
%  preln = tricoeffln(j1,j2,j3) + tricoeffln(j4,j5,j3) + tricoeffln(j1,j5,j6) + tricoeffln(j4,j2,j6);
%  termsln = preln/2 + facln(n+1) - ...
%    facln(alpha(1)-n) - facln(alpha(2)-n) - facln(alpha(3)-n) - ...
%    facln(n-beta(1)) - facln(n-beta(2)) - facln(n-beta(3)) -
%    facln(n-beta(4));
%
%  value = sum((-1).^n.*exp(termsln));
%if abs(value)<1e-10
%  value = 0;
%end
%elseif (ComputationMethod==1)
%  
%  n = max(beta):min(alpha);
%  d = tricoeff(a,b,c)*tricoeff(c,d,e)*tricoeff(a,e,f)*tricoeff(b,d,f);
%  t = factorial(alpha(1)-n).*factorial(alpha(2)-n).*factorial(alpha(3)-n).*...
%      factorial(n-beta(1)).*factorial(n-beta(2)).*factorial(n-beta(3)).*factorial(n-beta(4));
%  value = d*sum((-1).^n.*factorial(n+1)./t);
%
%if abs(value)<1e-10
%  value = 0;
%end
%elseif (ComputationMethod==2)

A = sort(alpha); A1 = A(1); A2 = A(2); A3 = A(3);
B = sort(beta); B1 = B(1); B2 = B(2); B3 = B(3); B4 = B(4);

pre = ...
  fac(A1-B2)*fac(A1-B3)*fac(A1-B4)/fac(A1-B1)/fac(B1+1) * ...
  fac(A2-B1)*fac(A2-B3)*fac(A2-B4)/fac(A2-B2)/fac(B2+1) * ...
  fac(A3-B1)*fac(A3-B2)*fac(A3-B4)/fac(A3-B3)/fac(B3+1);
pre = sqrt(pre*fac(B4+1));
sumfactor = (-1)^B4*bino(A1-B1,A1-B4)*bino(A2-B2,A2-B4)*bino(A3-B3,A3-B4);
pre = pre*sumfactor;

if (B4>=A1)
  value = pre;
else
  t = 1;
  s = t;
  for n = B4+1:A1
    t = -t*(n+1)/(n-B4)*(A1-n+1)/(n-B1)*...
           (A2-n+1)/(n-B2)*(A3-n+1)/(n-B3);
    s = s + t;
  end
  value = pre*s;
end
 
%end


%==========================================================
function v = bino(a,b)
v = prod((a-b+1:a)./(1:b));
return

function v = fac(a)
v = prod(1:a);
return

function ok = isint(x)
ok = x==floor(x);
return

function yes = istriangle(a,b,c)
yes = (a+b>=c) & (b+c>=a) & (c+a>=b);
return

%------------------------------------------------------------

function v = tricoeff(a,b,c)
v = sqrt(factorial(a+b-c)*factorial(a-b+c)/factorial(a+b+c+1)*factorial(-a+b+c));
return


function v = facln(x)
v = gammaln(x+1);
return


function v = tricoeffln(a,b,c)
v = (facln(a+b-c) + facln(a-b+c) + facln(-a+b+c) - facln(a+b+c+1));
return
