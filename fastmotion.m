% fastmotion  Computes linewidth parameters for fast-motion regime
%
%   lw = fastmotion(Sys,Field,tcorr)
%   [lw,mI] = fastmotion(...)
%
%   Computes FWHM Lorentzian linewidths
%   for the fast motional regime according to
%
%       FWHM = A + B*mI + C*mI^2
%
%   (one nucleus) and a similar formula for systems
%   with more than one nucleus. If present, the contribution
%   of the nuclear quadrupole is also included.
%
%   System    spin system structure
%      g      principal values of g matrix
%      gpa    g Euler angles (optional) [radians]
%      Nucs   nuclear isotope(s), e.g. '14N' or '63Cu'
%      A      principal values of A matrix [MHz]
%      Apa    A Euler angles (optional) [radians]
%      Q      principal values of Q tensor [MHz]
%   Field     center magnetic field [mT]
%   tcorr     isotropic rotational correlation time [seconds]
%
%   lw        all FWHM line widths [mT]
%   mI        mI quantum numbers for the lines, one line per row
%
%   Example:
%
%     System = struct('g',[2.0088 2.0064 2.0027],'Nucs','14N');
%     System.A = mt2mhz([7.59 5.95 31.76]/10);
%     [lw,mI] = fastmotion(System,350,1e-10)

% The undocumented syntax
%    [lw,mI,coeffs] = ...
% is for debugging purposes.

function [lw,mI,coeffs] = fastmotion(System,Field,tau20)

if (nargin==0), help(mfilename); return; end

if (nargin==0)
  % Di-t-butyl nitroxide
  System.Nucs = '14N';
  System.g = [2.0088 2.0064 2.0027];
  System.A = mt2mhz([7.59 5.95 31.76; 10 15 20]/10); % MHz
  Field = 350; % mT
  tau20 = 1e-10; % s
end

I = nucspin(System.Nucs);
nNucs = numel(I);

if numel(System.g)~=3
  error('Three values in System.g needed!');
end

if nNucs>0
  if numel(System.A)<nNucs*3
    error('Three values per nucleus in System.A needed!');
  end
end

if (tau20<1e-14)
  error('Correlation time too small!');
end

fullg = numel(System.g)==9;
if fullg
  g = System.g;
else
  if isfield(System,'gpa'), R = erot(System.gpa); else R = eye(3); end
  g = R*diag(System.g)*R.';
end

g0 = trace(g)/3;
g1 = reshape(g - eye(3)*g0,9,1);

if ~isfield(System,'Apa'); System.Apa = zeros(nNucs,3); end
if ~isfield(System,'Q'); System.Q = zeros(nNucs,3); end
if ~isfield(System,'Qpa'); System.Qpa = zeros(nNucs,3); end

for iNuc = 1:nNucs
  R = erot(System.Apa(iNuc,:));
  A = R*diag(System.A(iNuc,:))*R.';
  A0(iNuc) = trace(A)/3;
  A1(:,iNuc) = reshape(A - eye(3)*A0(iNuc),9,1);
end
if nNucs>0
  A1 = A1 * 1e6 *(2*pi); % angular frequency
end

% Coefficients A, B, C, and D
%-----------------------------------------------------
% Formulas taken from Atherton's book, p. 331, 332
% Scalar products of anistropic tensors:
%   Sum of squares of matrix elements.

cc = bmagn*Field*1e-3/(planck/2/pi); % g -> angular frequency
omega0 = g0*cc; % angular frequency

% Spectral densities j0 and j1
j0 = tau20;
j1 = tau20/(1+omega0^2*tau20^2);

% Anisotropy terms (all same units!)
if nNucs>0
  AA = A1.'*A1;
  gA = cc * g1.'*A1;
else
  AA = 0;
  gA = 0;
end
gg = cc^2 * g1.'*g1;
PP = sum(System.Q.^2,2).' * (1e6 *(2*pi))^2;

% g and A contributions
%------------------------------------------------------------
% Coefficients (in angular frequency units)
A = gg * (2/15*j0+1/10*j1) + sum(I.*(I+1).*diag(AA).') * (1/20*j0+7/60*j1);
B = gA * (4/15*j0+1/5*j1);
C = diag(AA).' * (1/12*j0-1/60*j1);
D = AA * (4/15*j0+1/10*j1);
D = D - tril(D);

A = convertcoeffs(A,g0);
B = convertcoeffs(B,g0);
C = convertcoeffs(C,g0);
D = convertcoeffs(D,g0);
coeffs.A = A;
coeffs.B = B;
coeffs.C = C;
coeffs.D = D;

% Nuclear quadrupole contribution
%-----------------------------------------------------------
% Formula taken from Atherton, p. 348. See also
%  - Freed/Fraenkel, JCP 39 (1963), 326-348, p. 338 (theory)
%  - Sillescu, Mol Phys 14 (1968), 381 (experiment, iodine I=5/2)
%  - Hudson, Luckhurst, Chem Rev 69 (1969) 191-225, p. 213 (review)
% qA: mI-independent, qB: mI^2, qC: mI^4
if nNucs>0
  II = I.*(I+1);
  qA = sum(j0/20 * PP.* II.*(II-1));
  qC = j0/20 * PP.* I.*(I+1).*2;
  qE = j0/20 * PP.* (-3);
  
  qA = convertcoeffs(qA,g0);
  qC = convertcoeffs(qC,g0);
  qE = convertcoeffs(qE,g0);
else
  qA = 0;
  qC = 0;
  qE = 0;
end
coeffs.qA = qA;
coeffs.qC = qC;
coeffs.qE = qE;

% Computation  of widths for all lines
%-----------------------------------------------------------
if nNucs>0
  mI = all_mI(I);
  for iLine = size(mI,1):-1:1
    mI_ = mI(iLine,:);
    lw(iLine) = A + sum(B.*mI_) + sum(C.*mI_.^2) + sum(sum(D.*(mI_'*mI_)));
    % add nqi contribution
    lw(iLine) = lw(iLine) + qA + sum(qC.*mI_.^2) + sum(qE.*mI_.^4);
  end
  lw = lw.';
else
  lw = A;
end

return
%============================================================
%============================================================
%============================================================

function c_mT = convertcoeffs(c_angfrq,g0)
%c = c_angfrq/1e6/(2*pi); % angular frequency -> MHz
c = c_angfrq/1e6/(1);
c = c/pi; % FWHM of Lorentzian is 1/pi/T2
c_mT = mhz2mt(c,g0); % -> mT
return


function mIall = all_mI(I)
nNucs = numel(I);

mIall = [];

nI = 2*I+1;
n = 1;
for iNuc = nNucs:-1:1
  mI = repmat(-I(iNuc):I(iNuc),n,1);
  n = n*nI(iNuc);
  mIall = [mI(:),repmat(mIall,nI(iNuc),1)];
end
return
