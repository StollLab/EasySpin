% convspec  Spectrum convolution with line shapes 
%
%   out = convspec(spec,df,fwhm)
%   out = convspec(spec,df,fwhm,deriv)
%   out = convspec(spec,df,fwhm,deriv,alpha)
%   out = convspec(spec,df,fwhm,deriv,alpha,phase)
%
%   Convolutes the spectral data array spec with a line
%   shape. df specifies the abscissa step, fwhm the
%   fwhm line width, deriv is the derivative (default 0),
%   alpha the shape parameter (1 is Gaussian, 0 is Loren-
%   tzian, default 1, between 0 and 1 a pseudovoigtian is
%   computed). phase is the phase of the Lorentzian component
%   (0 pure absorption, pi/2 pure dispersion).
%   out is the convoluted spectrum.
%
%   If the spectrum is more than 1D, different parameters
%   for each dimension can be defined. For 2D, e.g. df = 2
%   means that the abscissa step is 2 for both dimensions.
%   df = [2 3] means it is 2 for dim 1 and 3 for dim 2.
%   The parameters fwhm, deriv and alpha work similarly.
%   If fwhm is 0 for a certain dimension, convolution
%   along this dimension will be skipped.
%
%   Examples
%     spec = zeros(101,101);
%     spec(51,51) = 1; spec(30,60) = 3;
%     w = convspec(spec,.1,1);
%     pcolor(w);
%
%    convolutes in both dimensions with a Gaussian
%    with fwhm=1. The abscissa step is 0.1 for both
%    dimensions.
%
%     spec = zeros(100,100);
%     spec(51,51) = 1; spec(30,60) = 3;
%     w = convspec(spec,1,[3 6],[0 1],[0 1]);
%     surf(w);
%
%    convolutes in dim 1 with a Lorentzian fwhm=5 and
%    in dim 2 with a first derivative Gaussian fwhm=3.

function out = convspec(spec,steps,fwhm,deriv,alpha,phase)

if (nargin==0), help(mfilename); return; end

% Default values for optional input arguments
Default.alpha = 1; % pure Gaussian
Default.deriv = 0; % no derivative
Default.phase = 0; % pure absorption

switch nargin
  case 3
    deriv = Default.deriv;
    alpha = Default.alpha;
    phase = Default.phase;
  case 4
    alpha = Default.alpha;
    phase = Default.phase;
  case 5
    phase = Default.phase;
  case 6
  otherwise
    error('Wrong number of input arguments.');
end

if any(isnan(spec(:)))
  warning('The input data array (spec) contains NaN entries.');
end

if all(deriv~=-1) && all(deriv~=0) && all(deriv~=1) && all(deriv~=2)
  error('deriv must be 0, 1 or 2.');
end

if any(alpha<0) || any(alpha>1)
  error('alpha must satisfy 0 <= alpha <= 1.');
end

% Remove all singleton dimensions. The original array shape is restored at the end.
sz_original = size(spec);
sz = sz_original;
sz(sz==1) = [];
if numel(sz)==1
  spec = reshape(spec,[sz,1]);
else
  spec = reshape(spec,sz);
end
nDims = length(sz);

if (nDims>2) && (numel(spec)>200^3)
  warning('Very big data array! convspec might take very long!');
end

expandfull = @(in,n_)[in(:).' in(ones(1,n_-length(in)))];

steps = expandfull(steps,nDims);
fwhm = expandfull(fwhm,nDims);
deriv = expandfull(deriv,nDims);
alpha = expandfull(alpha,nDims);
phase = expandfull(phase,nDims);

NN = 2*sz + 1;       % extended range for each dimension
mid = round(NN/2)+1; % range midpoint, for line center
fwhm = fwhm./steps;  % line width in terms of x-axis increments
steps(fwhm==0) = 1;

% Error if linewidth for first-derivative convolution is too small
if any(fwhm<1 & fwhm>0 & deriv==1)
  error(['The linewidth is smaller (%gx) than the increment. ' ...
    'The convolution with such a narrow derivative lineshape does not work. ' ...
    'Either increase linewidth or increase axis resolution.'],min(fwhm(fwhm~=0)));
end

% Determine ifft of line shape for each dimension
for i = nDims:-1:1
  Range{i} = 1:sz(i);
  if fwhm(i)>0
    Line = lshape(1:NN(i),mid(i),fwhm(i),deriv(i),alpha(i),phase(i)).';
    LineDecay{i} = ifft(Line([mid(i):end 1:mid(i)-1]));
  else
    LineDecay{i} = ones(1,NN(i))/NN(i);
  end
end

% Convolve with general outer product of line decays
out = fftn(ifftn(spec,[NN 1]) .* xouter(LineDecay{:}));
out = prod(NN)/prod(steps.^deriv) * out;

% Extract original range and reshape to original shape
out = reshape(out(Range{:}),sz_original);

% If the input was real, the output has to be, too.
% Remove imaginary numerical noise.
if isreal(spec), out = real(out); end

return
%=======================================================================


%=======================================================================
% xouter  general outer product
%
%   A = xouter(x1,x2,...)
%
%   Computes the outer product of the vectors
%   x1, x2, etc. A has as many dimensions as
%   there are input vectors. E.g. three inputs;
%
%   A(i,j,k) = x1(i)*x2(j)*x3(k)

function out = xouter(varargin)

switch nargin
case 0, out = zeros(0); return;
case 1, out = varargin{1}(:); return;
end

nDims = length(varargin);

for i=length(varargin):-1:1
  n(i) = numel(varargin{i});
end

out = ones(n);

for i = 1:nDims
  x = varargin{i}(:); % Extract and reshape as a vector.
  lengths = n;
  lengths(i) = []; % Remove i-th dimension
  x = reshape(x(:,ones(1,prod(lengths))),[n(i) lengths]); % Expand x
  out = out .* permute(x,[2:i 1 i+1:nDims]); % Permute to i'th dimension
end

return
%=======================================================================
