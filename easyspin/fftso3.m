% fftso3   Fast Fourier transform over the rotation group SO(3)
%
%   c = fftso3(fcn,Lmax)
%   [LMK,vals] = fftso3(fcn,Lmax)
%   [c,LMK,vals] = fftso3(fcn,Lmax)
%   ___ = fftso3(fcn,Lmax,mode)
%
% Takes the function fcn(alpha,beta,gamma) defined over the Euler angles
% 0<=alpha<2*pi, 0<=beta<=pi, 0<=gamma<2*pi and calculates the coefficients
% c^L_MK of its tuncated expansion in terms of Wigner D-functions, D^L_MK:
%
%    fcn(alpha,beta,gamma) = sum_(L,M,K) c^L_MK * D^L_MK(alpha,beta,gamma)
%
% where L = 0,..,Lmax, M = -L,..,L, K = -L,..,L.
%
% Input:
%   fcn   function handle for function fcn(alpha,beta,gamma), with angles
%         in units of radians; it should be vectorized and accept array
%         inputs for alpha and gamma
%   Lmax  maximum L for expansion
%   mode  beta grid mode, either 'GL' for Gauss-Legendre (Lmax+1 points)(default)
%         or 'DH' for Driscoll-Healy (2*Lmax+2 points)
%
% Output:
%   c     (Lmax+1)-element cell array of Fourier coefficients, L = 0 .. Lmax
%         c{L+1} contains the (2L+1)x(2L+1)-dimensional matrix of expansion
%         coefficients c^L_MK, with M and K ordered in ascending order: M,K =
%         -L,-L+1,..,L.
%   LMK   L, M, and K values of the non-zero coefficients
%   vals  values of the non-zero coefficients

% References:
% 1) D. K. Maslen, D. N. Rockmore
%    Generalized FFT's - A survey of some recent results
%    Groups and Computation II
%    DIMACS Series in Discrete Mathematics and Theoretical Computer
%    Science, vol.28
%    L. Finkelstein, W. M. Kantor (editors)
%    American Mathematical Society, 1997
%    https://doi.org/10.1007/s00041-008-9013-5
% 2) P. J. Kostelec, D. N. Rockmore
%    FFTs on the Rotation Group
%    J.Fourier.Anal.Appl. 2008, 14, 145/179
%    https://doi.org/10.1007/s00041-008-9013-5
% 3) Z. Khalid, S. Durrani, R. A. Kennedy, Y. Wiaux, J. D. McEwen
%    Gauss-Legendre Sampling on the Rotation Group
%    IEEE Signal Process. Lett. 2016, 23, 207-211
%    https://doi.org/10.1109/LSP.2015.2503295

function varargout = fftso3(fcn,Lmax,betamode)

if nargin==0
  help(mfilename);
end

if nargin~=2 && nargin~=3
  error('Two inputs are required (fcn, Lmax).');
end

if nargin<3
  betamode = 'GL';
end

if numel(Lmax)~=1 || mod(Lmax,1) || Lmax<0
  error('Second input (Lmax) must be a non-negative integer.');
end

B = Lmax + 1; % bandwidth (0 <=L < B)

% Set up sampling grid over (alpha,gamma)
alpha = 2*pi*(0:2*B-2)/(2*B-1); % radians
gamma = 2*pi*(0:2*B-2)/(2*B-1); % radians

% Get knots and weights for integration over beta
switch betamode
  case 'DH' % Driscoll-Healy
    beta = pi*((0:2*B-1)+1/2)/(2*B); % radians
    nbeta = numel(beta);
    % Calculate beta weights
    q = 0:B-1; % summation index
    w = zeros(1,nbeta);
    for b = 1:nbeta
      w(b) = 2/B * sin(beta(b)) * sum(sin((2*q+1)*beta(b))./(2*q+1));
    end
  case 'GL' % Gauss-Legendre
    [z,w] = gauleg(-1,1,B);
    beta = acos(z);
end
nbeta = numel(beta);

% Calculate 2D inverse FFT along alpha and gamma, one for each beta
[alpha,gamma] = ndgrid(alpha,gamma); % 2D grid for function evaluations
S = cell(1,nbeta);
for b = 1:nbeta
  f = fcn(alpha,beta(b),gamma); % evaluate function over (alpha,gamma) grid
  S_ = (2*pi)^2*ifft2(f);
  S{b} = fftshift(S_); % reorder to M,K = -(B-1), ..., 0, ..., B-1
end

% Calculate Wigner d-function values for all L and beta, with matrix exponential
d = cell(Lmax+1,nbeta);
for L = 0:Lmax
  v = sqrt((1:2*L).*(2*L:-1:1))/2i; % off-diagonals of Jy matrix
  Jy = diag(v,-1) - diag(v,1); % Jy matrix (M,K = -Lmax,-Lmax+1,...,+Lmax)
  [U,mj] = eig(Jy);
  mj = diag(mj).'; % equal to -L:L
  for b = 1:nbeta
    e = exp(-1i*beta(b)*mj);
    d{L+1,b} = (U .* e) * U'; % equivalent to U*diag(e)*U', but slightly faster
  end
end

% Calculate Fourier coefficients
c = cell(Lmax+1,1);
for L = 0:Lmax
  c_ = zeros(2*L+1);
  idx = (Lmax+1)+(-L:L); % to access the range M,K = -L:L in S
  for b = 1:nbeta % run over all beta values
    c_ = c_ + w(b)*d{L+1,b}.*S{b}(idx,idx);
  end
  c{L+1} = (2*L+1)/(8*pi^2)*c_;
end

% Remove numerical noise
maxcoeff = cellfun(@(x)max(abs(x),[],'all'),c);
threshold = 1e-13*max(maxcoeff);
for L = 0:Lmax
  idx = abs(c{L+1})<threshold;
  c{L+1}(idx) = 0;
end

% Return full matrices or list of non-zero coefficients
returnList = nargout>1;
if returnList
  % Convert to [L M K value] list
  LMK = [];
  vals = [];
  for L = 0:Lmax
    [idxK,idxM,vals_] = find(c{L+1}.'); % transpose to assure lexicographic order of L,M,K
    nNonzero = numel(vals_);
    if nNonzero>0
      LMK = [LMK; repmat(L,nNonzero,1) idxM-L-1 idxK-L-1];
      vals = [vals; vals_];
    end
  end
end

% Organize output
switch nargout
  case 1, varargout = {c};
  case 2, varargout = {LMK,vals};
  case 3, varargout = {c,LMK,vals};
end

end

% Calculate Gauss-Legendre grid points and weights over inteval [x1,x2]
% Press et al, Numerical Recipes in C, 2nd ed., p.152
function [x,w] = gauleg(x1,x2,n)

m = (n+1)/2;
xm = (x2+x1)/2;
xl = (x2-x1)/2;
for i = 1:m
  z = cos(pi*(i-0.25)/(n+0.5));
  z1 = inf;
  while abs(z-z1)>eps
    p1 = 1;
    p2 = 0;
    for j = 1:n
      p3 = p2;
      p2 = p1;
      p1 = ((2*j-1)*z*p2 - (j-1)*p3)/j;
    end
    pp = n*(z*p1-p2)/(z^2-1);
    z1 = z;
    z = z1 - p1/pp;
  end
  x(i) = xm - xl*z;
  x(n+1-i)= xm + xl*z;
  w(i) = 2*xl/((1-z^2)*pp^2);
  w(n+1-i) = w(i);
end

end
