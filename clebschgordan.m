% clebschgordan  Clebsch-Gordan coefficient
%
%   v = clebschgordan(j1,j2,j,m1,m2,m)
%   v = clebschgordan(jm1,jm2,jm)
%   v = clebschgordan(j,m)
%
%   Returns the Clebsch-Gordan coefficient
%
%      (j1,j2,m1,m2|j1,j2,j,m)
%
%   Definitons for alternative input forms
%   a)  jm1 = [j1 m1] jm2 = [j2 m2], jm = [j m]
%   b)  j = [j1 j2 j], m = [m1 m2 m]

function value = clebschgordan(varargin)

if (nargin==0), help(mfilename); return; end

% Parse input
switch (nargin)
  case 6
    [j1,j2,j,m1,m2,m] = deal(varargin{:});
  case 3
    [jm1,jm2,jm] = deal(varargin{:});
    j1 = jm1(1); m1 = jm1(2);
    j2 = jm2(1); m2 = jm2(2);
    j = jm(1); m = jm(2);
  case 2
    [j,m] = deal(varargin{:});
    j1 = j(1); j2 = j(2); j = j(3);
    m1 = m(1); m2 = m(2); m = m(3);
  otherwise
    error('Wrong number of parameters!');
end

% Compute coefficient via 3-j symbol
value = (-1)^(j1-j2+m) * sqrt(2*j+1) * wigner3j([j1,j2,j],[m1,m2,-m]);

return
