% voigtian  Voigtian line shape function 
%
%    y = voigtian(x,x0,fwhm);
%    y = voigtian(x,x0,fwhm,deriv);
%
%    Computes the convolution of a Gaussian and a
%    Lorentzian line shape function.
%
%    x     abscissa vector
%    x0    line shape centre
%    fwhm  [fwhm_Gauss fwhm_Lorentz]
%    deriv 0 (absorption), 1 (first), 2 (second
%          derivative)
%
%    y     Voigtian line shape, normalized to
%          integral 1

function y = voigtian(x,x0,fwhm,deriv)

if nargin==0, help(mfilename); return; end

if (nargin<3) || (nargin>4), error('Wrong number of input arguments!'); end

if numel(fwhm)~=2
  error('fwhm must contain 2 values!');
end

if nargin==3, deriv=0; end

derivG = 0;
derivL = 0;
if fwhm(1)>fwhm(2)
  derivG = deriv;
else
  derivL = deriv;
end

% Compute Gaussian and Lorentzian
% one of the two lineshapes must be centered in range to yield
% correct final lineshape position
y1 = gaussian(x,x0,fwhm(1),derivG);
y2 = lorentzian(x,x(ceil(numel(x)/2)),fwhm(2),derivL);

% Convolution
y = conv(y1,y2);

% Take central part
m = ceil(length(y1)/2);
y = y(m:end-length(y1)+m);

% Re-normalize
dx = x(2)-x(1);
y = y*dx;

return

clear
x0 = -5; ha = 0;
fwhm = [2 1]/2;
x = linspace(-10,10,5001);
y = voigtian(x,x0,fwhm,ha);
y1 = gaussian(x,x0,fwhm(1),ha);
y2 = lorentzian(x,x0,fwhm(2),ha);

plot(x,y,'r',x,y1,'b',x,y2,'g');
