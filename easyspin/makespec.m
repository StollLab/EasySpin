% makespec   Construct stick spectrum from peak data 
%
%  spec = makespec(Range,nPoints,Pos)
%  spec = makespec(Range,nPoints,Pos,Amp)
%  [x,spec] = makespec(...)
%
%  Constructs a stick spectrum from peak positions and intensities.
%
%   Range    limits of the abscissa [minX maxX]
%   nPoints  length of spectral vector
%   Pos      array of peak positions
%   Amp      array of peak amplitudes
%            if omitted all amplitudes are set to 1
%
%   spec     spectrum
%   x        abscissa vector
%
%   Pos and Amp must contain the same number of elements

function varargout = makespec(Range,nPoints,Positions,Amplitudes)

if (nargin==0), help(mfilename); return; end

if (nargin==3)
  Amplitudes = ones(numel(Positions),1);
end

if (numel(Amplitudes)~=numel(Positions))
  error('The position and amplitude vectors must have the same number of elements!');
end

if (numel(Range)~=2) | ~isreal(Range)
  error('Range must be a 2-element array [minX maxX].');
end
   
if (Range(1)>=Range(2))
  error('Range(2) must be larger than Range(1).');
end

idxPositions = (nPoints-1) * (Positions-Range(1))/diff(Range);
idxPositions = 1 + fix(idxPositions);
OutOfRange = (idxPositions<1) | (idxPositions>nPoints);

if any(OutOfRange)
  warning('Some peaks fall outside the specified range.');
  idxPositions(OutOfRange) = [];
  Amplitudes(OutOfRange) = [];
end

spec = full(sparse(1,idxPositions,Amplitudes,1,nPoints));

switch nargout
case {0,1},
  varargout = {spec};
case 2,
  x = linspace(Range(1),Range(2),nPoints);
  varargout = {x,spec};
end

return

%---------------------------------------------------
function makespec_test

clear, close all
N = 20;
pos = rand(1,N); amp = rand(1,N);
[x,s] = makespec([-0.2,1.2],2000,pos,amp);
subplot(2,1,1); plot(x,s)

N = 200e3;
pos = randn(1,N); amp = rand(1,N);
[x,s] = makespec([-1,2],2000,pos,amp);
subplot(2,1,2); plot(x,s)

return
