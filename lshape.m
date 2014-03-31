% lshape  Gaussian and Lorentzian line shapes 
%
%   y = lshape(x,x0,fwhm)
%   y = lshape(x,x0,fwhm,diff)
%   y = lshape(x,x0,fwhm,diff,alpha)
%   y = lshape(x,x0,fwhm,diff,alpha,phase)
%
%   General normalized line shape function, a linear
%   combination of Lorentzian and Gaussian lineshapes.
%
%   Input:
%   - x: Abscissa vector
%   - x0: Center of the linshape function
%   - fwhm: Full width at half height. If alpha is
%     given, and Gaussian and Lorentzian components should
%     have different fwhm, specify a vector
%         [fwhmGauss fwhmLorentz]
%   - diff: Derivative. 0 is no derivative, 1 first,
%     2 second, -1 the integral with -infinity as lower
%     limit. If omitted, 0 is the default.
%   - alpha: Weighting factor for linear combination
%       alpha*Gaussian + (1-alpha)*Lorentzian.
%     Gaussian is default.
%   - phase: phase (0 pure absorption, pi/2 pure dispersion)
%
%   Output:
%   - y: Vector of function values for abscissa x

function y = lshape(x,x0,fwhm,varargin)

if (nargin==0), help(mfilename); return; end

% Functions are normalized for diff=0: the integral
% from -inf to inf is 1. Derivatives of these
% functions are with respect to x.

% parse input parameters
diff = 0; % default
alphaGauss = 1; % default 
phase = 0; % default
switch nargin
 case 4
  diff = varargin{1};
 case 5
  diff = varargin{1};
  alphaGauss = varargin{2};
 case 6
  diff = varargin{1};
  alphaGauss = varargin{2};
  phase = varargin{3};
end

if any(fwhm<=0) || any(~isreal(fwhm))
  error('fwhm must be positive and real!');
end

if numel(fwhm)>2
  error('fwhm must contain 1 or 2 elements!');
end

if numel(fwhm)==2
  fwhmG = fwhm(1);
  fwhmL = fwhm(2);
elseif numel(fwhm)==1
  fwhmG = fwhm(1);
  fwhmL = fwhm(1);
end

% compute Gaussian component, if demanded
if (alphaGauss~=0)
  [GaussAbs,GaussDisp] = gaussian(x,x0,fwhmG,diff);
  if (phase~=0)
    Gauss = GaussAbs*cos(phase) +  GaussDisp*sin(phase);
  else
    Gauss = GaussAbs;
  end
else
  Gauss = 0;
end

% compute Lorentzian component, if demanded
if (alphaGauss~=1)
  [LorentzAbs,LorentzDisp] = lorentzian(x,x0,fwhmL,diff);
  if (phase~=0)
    Lorentz = LorentzAbs*cos(phase) +  LorentzDisp*sin(phase);
  else
    Lorentz = LorentzAbs;
  end
else
  Lorentz = 0;
end

% build hybrid line shape, also known as pseudo-Voigt function
y = alphaGauss*Gauss + (1-alphaGauss)*Lorentz;

return
