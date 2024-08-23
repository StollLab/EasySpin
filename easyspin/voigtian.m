% voigtian  Voigtian line shape function 
%
%   ya = voigtian(x,x0,fwhmGL)
%   ya = voigtian(x,x0,fwhmGL,diff)
%   ya = voigtian(x,x0,fwhmGL,diff,phase)
%   [ya,yd] = voigtian(...)
%
%    Computes an area-normalized Voigt line shape function, which is
%    the convolution of a Gaussian and a Lorentzian line shape function.
%
%   Input:
%   - x:     abscissa vector
%   - x0:    center of the lineshape function
%   - fwhmGL:  full widths at half height [fwhm_Gauss fwhm_Lorentz]
%   - diff:  derivative. 0 is no derivative, 1 first,
%            2 second and so on, -1 the integral with -infinity
%            as lower limit. 0 is the default.
%   - phase: phase rotation, mixes absorption and dispersion.
%            phase=pi/2 puts dispersion signal into ya
%
%   Output:
%   - ya:   absorption function values for abscissa x
%   - yd:   dispersion function values for abscissa x

function [y1,y2] = voigtian(x,x0,fwhmGL,diff,phase)

if nargin==0, help(mfilename); return; end

if nargin<3, error('Not enough input arguments!'); end
if nargin>5, error('Too many input arguments!'); end

if ~isnumeric(x0) || numel(x0)~=1
  error('x0 (2nd input) must be a single number.');
end

if ~isnumeric(fwhmGL) || numel(fwhmGL)~=2
  error('fwhmGL (3rd input) must contain two numbers!');
end

if ~isreal(fwhmGL) || any(fwhmGL<0)
  error('Both numbers in fwhmGL (3rd input) must be non-negative.');
end

if all(fwhmGL==0)
  error('At least one value in fwhmGL (3rd input) must be nonzero.');
end

if nargin<4, diff = 0; end
if nargin<5, phase = 0; end

calculateQuadrature = nargout>1;

fwhmG = fwhmGL(1);
fwhmL = fwhmGL(2);

% For the convolution, one of the two lineshapes must be centered in range to
% yield the correct final lineshape position
x0center = x(ceil(numel(x)/2));

% Use the wider lineshape to include position, phase, and derivative.
% Use the narrower one for convolution only.
if fwhmG>fwhmL
  diffG = diff;
  phaseG = phase;
  x0G = x0;
  diffL = 0;
  phaseL = 0;
  x0L = x0center;  
else
  diffG = 0;
  phaseG = 0;
  x0G = x0center;  
  diffL = diff;
  phaseL = phase;
  x0L = x0;
end

% Compute Gaussian and Lorentzian
if calculateQuadrature
  if fwhmG>0
    [yG,yGquad] = gaussian(x,x0G,fwhmG,diffG,phaseG);
  end
  if fwhmL>0
    [yL,yLquad] = lorentzian(x,x0L,fwhmL,diffL,phaseL);
  end
else
  if fwhmG>0
    yG = gaussian(x,x0G,fwhmG,diffG,phaseG);
    yGquad = [];
  end
  if fwhmL>0
    yL = lorentzian(x,x0L,fwhmL,diffL,phaseL);
    yLquad = [];
  end
end

doConvolution = fwhmL>0 && fwhmG>0;
if doConvolution
  
  % Convolution
  y1 = conv(yG,yL,'same');
  if calculateQuadrature
    if fwhmL>fwhmG
      y2 = conv(yG,yLquad,'same');
    else
      y2 = conv(yGquad,yL,'same');
    end
  end
    
  % Re-normalize
  dx = x(2)-x(1);
  y1 = y1*dx;
  if calculateQuadrature
    y2 = y2*dx;
  end
  
else
  
  % Skip convolution in the case of a pure Gaussian or pure Lorentzian
  if fwhmL>0
    y1 = yL;
    y2 = yLquad;
  else
    y1 = yG;
    y2 = yGquad;
  end
  
end

end
