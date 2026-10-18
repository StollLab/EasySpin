% makespec   Construct stick spectrum from peak data
%
%  spec = makespec(Range,nPoints,Pos)
%  spec = makespec(Range,nPoints,Pos,Amp)
%  spec = makespec(Range,nPoints,Pos,Amp,Method)
%  [x,spec] = makespec(...)
%  [x,spec,nOut] = makespec(...)
%
%  Constructs a stick spectrum from peak positions and intensities.
%
%   Range    limits of the abscissa [minX maxX]
%   nPoints  length of spectral vector
%   Pos      array of peak positions
%   Amp      array of peak amplitudes, or a single amplitude for all peaks
%            if omitted or empty, all amplitudes are set to 1
%   Method   how peaks are assigned to points of the abscissa
%            'round'  nearest point (default)
%            'floor'  nearest point at or below the peak position
%            'split'  amplitude split between the two neighboring points,
%                     in proportion to their proximity to the peak
%
%   spec     spectrum
%   x        abscissa vector
%   nOut     number of peaks outside the range, [nBelow nAbove]
%            if requested, no warning is issued for out-of-range peaks
%
%   Peaks with NaN positions are ignored.

function varargout = makespec(Range,nPoints,Positions,Amplitudes,Method)

if nargin==0, help(mfilename); return; end

if nargin<4 || isempty(Amplitudes)
  Amplitudes = 1;
end

if nargin<5 || isempty(Method)
  Method = 'round';
end

if ~isscalar(Amplitudes) && numel(Amplitudes)~=numel(Positions)
  error('The position and amplitude vectors must have the same number of elements!');
end

if (numel(Range)~=2) || ~isreal(Range)
  error('Range must be a 2-element array [minX maxX].');
end

if Range(1)>=Range(2)
  error('Range(2) must be larger than Range(1).');
end

Methods = {'round','floor','split'};
if ~(ischar(Method) || isstring(Method)) || ~any(strcmpi(Method,Methods))
  error('Method must be one of ''round'', ''floor'', or ''split''.');
end
Method = lower(char(Method));

% Convert positions to zero-based fractional indices into the spectral vector
u = (Positions(:)-Range(1))/diff(Range)*(nPoints-1);
Amplitudes = Amplitudes(:).*ones(size(u));

switch Method
  case 'round'
    idx = 1 + round(u);
    below = idx<1;
    above = idx>nPoints;
  case 'floor'
    idx = 1 + floor(u);
    below = idx<1;
    above = idx>nPoints;
  case 'split'
    below = u<0;
    above = u>nPoints-1;
end

keep = ~below & ~above & ~isnan(u);
nOut = [nnz(below) nnz(above)];

if strcmp(Method,'split')
  % Distribute each amplitude linearly over the two neighboring points
  u = u(keep);
  a = Amplitudes(keep);
  i0 = floor(u);
  f = u - i0;
  atEnd = i0==nPoints-1;  % peak exactly at Range(2)
  i0(atEnd) = nPoints-2;
  f(atEnd) = 1;
  idx = [i0+1; i0+2];
  Amplitudes = [(1-f).*a; f.*a];
else
  idx = idx(keep);
  Amplitudes = Amplitudes(keep);
end

if any(nOut) && nargout<3
  warning('Some peaks fall outside the specified range.');
end

spec = accumarray(idx(:),Amplitudes(:),[nPoints 1]).';

switch nargout
  case {0,1}
    varargout = {spec};
  case {2,3}
    x = linspace(Range(1),Range(2),nPoints);
    varargout = {x,spec,nOut};
end

end
