% fftso3   Fast Fourier transform over the rotation group SO(3)
%
%   c = fftso3(fcn,Lmax)
%   c = fftso3(fcn,Lmax,nrm)
%
% Takes the function fcn(aplha,beta,gamma) defined over the Euler angles
% 0<=alpha<2*pi, 0<=beta<=pi, 0<=gamma<2*pi and calculates the coefficients c^L_MK
% of its truncated Fourier expansion in terms of Wigner D-functions D^L_MK:
%
%    fcn(alpha,beta,gamma) = sum_(L,M,K) f^L_MK * D^L_MK(alpha,beta,gamma)
%
% where L = 0..Lmax, M = -L..L, K = -L..L.
%
% Input:
%   fcn   function handle for function fcn(alpha,beta,gamma), with angles
%         in units of radians
%   Lmax  maximum L for Wigner expansion
%   nrm   true: use normalized Wigner functions sqrt((2L+1)/8pi^2)*D^L_MK
%         false: use standard un-normalized Wigner functions D^L_MK
%           (default)
%
% Output:
%   c     (Lmax+1)-element cell array of Fourier coefficients, L = 0 .. Lmax
%         c{L+1} contains the (2L+1)x(2L+1) matrix of Fourier coefficients
%         f^L_MK, with M and K ordered in ascending order: M,K =
%         -L,-L+1,..,L

% References:
% 1) P. J. Kostelec, D. N. Rockmore
%    FFTs on the Rotation Group
%    J.Fourier.Anal.Appl. 2008, 14, 145/179
%    https://doi.org/10.1007/s00041-008-9013-5
% 2) D.-M. Lux, C. Wülker, G. S. Chirikjian,
%    Parallelization of the FFT on SO(3)
%    arXiv:1808.00896 [cs.DC]
%    https://arxiv.org/abs/1808.00896

function c = fftso3(fcn,Lmax,nrm)

if nargin==0
  help(mfilename);
end

if nargin<3
  nrm = false;
end

B = Lmax+1; % bandwidth (0 <=L < B)

% Set up sampling grid over (alpha,beta,gamma)
k = 0:2*B-1;
alpha = 2*pi*k/(2*B); % radians
beta = pi*(2*k+1)/(4*B); % radians
gamma = 2*pi*k/(2*B); % radians

% Calculate 2D inverse FFT along alpha and gamma, one for each beta
S = cell(1,2*B);
for j = 0:2*B-1
  f = fcn(alpha.',beta(j+1),gamma); % evaluate function over (alpha,gamma) grid
  S_ = (2*B)^2*ifft2(f); % remove the prefactor included by ifft2()
  S_ = fftshift(S_); % reorder to M1,M2 = -B, -B+1, ..., 0, ..., B-1
  S{j+1} = S_(2:end,2:end); % drop first row and column (M1,M2 = -B)
end

% Calculate beta weights
i = 0:B-1; % summation index
wB = zeros(1,2*B);
for j = 0:2*B-1
  wB(j+1) = 2*pi/B^2 * sin(beta(j+1)) * sum(sin((2*i+1)*beta(j+1))./(2*i+1));
end

% Calculate all Wigner d-function values for all L and beta, using matrix exponential
d = cell(Lmax+1,2*B);
for L = 0:Lmax
  v = sqrt((1:2*L).*(2*L:-1:1))/2; % off-diagonals of Jy matrix (without 1/i)
  Jy = diag(v,-1) - diag(v,1); % assemble Jy matrix (M1,M2 = -Lmax,-Lmax+1,...,Lmax)
  for j = 0:2*B-1
    d{L+1,j+1} = expm(-beta(j+1)*Jy); % (i dropped since 1/i dropped in Jy)
  end
end

% Calculate Fourier coefficients
c = cell(1,Lmax+1);
for L = 0:Lmax
  c_ = zeros(2*L+1);
  idx = (Lmax+1)+(-L:L);
  for j = 0:2*B-1 % run over all beta values
    c_ = c_ + wB(j+1)*d{L+1,j+1}.*S{j+1}(idx,idx);
  end
  c{L+1} = (2*L+1)/(8*pi*B)*c_;
end

% Normalize
if nrm
  for L = 0:Lmax
    N = sqrt((2*L+1)/(8*pi^2));
    c{L+1} = c{L+1}/N;
  end
end

% Zero coefficients that are below threshold
maxcoeff = cellfun(@(x)max(abs(x),[],'all'),c);
threshold = 1e-13*max(maxcoeff);
for L = 0:Lmax
  idx = abs(c{L+1})<threshold;
  c{L+1}(idx) = 0;
end

end
