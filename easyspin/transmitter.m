% transmitter      Simulation of/compensation for the effect of transmitter
%                  nonlinearity on a pulse
%
%  signal = transmitter(signal0,Ain,Aout,'simulate')
%  signal = transmitter(signal0,Ain,Aout,'compensate')
%
%  If the option 'simulate' is selected, transmitter() simulates the effect
%  of transmitter nonlinearity provided as output amplitude (Aout) as a
%  function of input amplitude (Ain) on the input signal signal0.
%
%  If the option 'compensate' is selected, transmitter() adapts the signal
%  provided in signal0 to compensate for the defined transmitter nonlinearity.
%
%  Input:
%   - signal0   = input signal vector, for 'simulate' the amplitude of signal0
%                 refers to the amplitude scale provided in Ain, for 'compensate'
%                 it refers to the amplitude scale Aout
%   - Ain       = input amplitudes
%   - Aout      = output amplitudes
%   - 'simulate'/'compensate'
%
%  Output:
%   - signal    = signal with simulated amplitude compression or signal compensated
%                 for the transmitter nonlinearity
%

function signal = transmitter(signal0,Ain,Aout,option,varargin)

% Input argument check
% ----------------------------------------------------------------------- %
if nargin==0
  help(mfilename);
  return
end

if nargin~=4
  error('transmitter() requires 4 input arguments.')
end

% Check signal dimensions
dim = size(signal0);
if min(dim)>1
  error('Input must be a row or column vector.');
end
if numel(dim)>2
  error('Input must be a row or column vector.');
end

% Check input and output amplitudes
if numel(Ain)~=numel(Aout)
  error('The input and output amplitude vectors need to have the same length.')
end

if strcmp(option,'simulate') && max(abs(signal0))>max(Ain)
  error('The input signal amplitude needs to be defined on the provided input amplitude scale of the transmitter.')
end
if strcmp(option,'compensate') && max(abs(signal0))>max(Aout)
  error('The amplitude of the desired output signal needs to be within the provided output amplitude scale of the transmitter.')
end

% Default options
if nargin==5
  Opt = varargin{1};
else
  Opt = struct();
end
if ~isfield(Opt,'N') || isempty(Opt.N)
  Opt.N = 4;
end

% Transmitter power transfer curve
% ----------------------------------------------------------------------- %
Ain = Ain(:);
Aout = Aout(:);

% Polynomial fit (constrained to origin)
% Vandermonde matrix
V(:,Opt.N-1) = Ain;
for j = Opt.N-2:-1:1
  V(:,j) = Ain.*V(:,j+1);
end
coeff = [V\Aout; 0];
Aout = polyval(coeff,Ain);

% Transmitter nonlinearity simulation and compensation
% ----------------------------------------------------------------------- %
if verLessThan('Matlab','7.13') % R2011b
  % because of calls to griddedInterpolant
  error('Your Matlab version does not support this feature. At least R2011b is needed.');
end

switch option
  case 'simulate'
    
    F = griddedInterpolant(Ain,Aout,'spline');
    
  case 'compensate'
    
    [OutputAmplitude,uniqueind] = unique(Aout);
    InputAmplitude = Ain(uniqueind);
    F = griddedInterpolant(OutputAmplitude,InputAmplitude,'spline');
    
  otherwise
    
    error('Requested calculation option not available. Available options are ''simulate'' and ''compensate''.');
    
end

% Compressed/compensated signal
if isreal(signal0)
  signal = sign(signal0).*F(abs(signal0));
else
  signal = sign(real(signal0)).*F(abs(real(signal0))) + ...
    1i*sign(imag(signal0)).*F(abs(imag(signal0)));
end

return