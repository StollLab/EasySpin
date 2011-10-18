% fdaxis  Frequency domain axis 
%
%   xf = fdaxis(dT,N)
%   xf = fdaxis(xt)
%
%   Returns a vector xf containing the frequency-
%   domain axis of the FFT of a N-point time-
%   domain vector xt sampled with period dT.
%   You can either specify xt or dT and N.
%
%   xf has 0 (DC component) in the center. Units:
%   [xf] = 1/[dT]. Time yields frequency, NOT
%   angular frequency.

function FreqAxis = fdaxis(varargin)

if (nargin==0), help(mfilename); return; end

switch nargin
case 1
  TimeAxis = varargin{1}(:);
  dT = TimeAxis(2) - TimeAxis(1);
  N = length(TimeAxis);
case 2
  dT = varargin{1};
  N = varargin{2};
otherwise
  error('Wrong number of input arguments!');
end

% dT: sampling time in time domain
% N:  number of samples in time domain

NyquistFrequency = 1/(2*dT);
UnitAxis = 2/N * ((0:N-1)-fix(N/2));
FreqAxis = NyquistFrequency * UnitAxis;
% DC component is in the center.

return
