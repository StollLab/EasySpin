% addnoise   Add noise to a signal
%
%   yn = addnoise(y,SNR,'f')
%   yn = addnoise(y,SNR,'n')
%   yn = addnoise(y,SNR,'u')
%
%   Adds noise to y to give yn with signal-to-noise ratio SNR. The type of
%   noise is specified as third input:
%   'f': 1/f noise
%   'n': normal (Gaussian) distribution
%   'u': uniform distribution
%
%   The SNR is defined as the ratio of signal amplitude to the
%   standard deviation of the noise amplitude distribution.
%
%   Example:
%     x = linspace(-1,1,1001);
%     y = lorentzian(x,0,0.3,1);
%     yn = addnoise(y,20,'n');
%     plot(x,yn);

function yn = addnoise(y,SNR,noisemodel)

if nargin==0, help(mfilename); return; end

if nargin<2
  error('Second input must be the desired signal-to-noise ratio (SNR).');
end
if nargin<3
  error('Third input must be the desired noise type (''f'', ''n'', or ''u'').');
end

if numel(SNR)~=1 || SNR<=0
  error('The second input, the signal-to-noise ratio (SNR), must be a positive number.');
end

if ~any(noisemodel=='unf')
  error('The third input, the noise model, must be ''u'', ''n'', or ''f''.');
end

dims = size(y);

complexData = ~isreal(y); % check for complex data

switch noisemodel
  case 'u' % uniform distribution [-0.5,0.5]
    rnoise = rand(dims)-0.5;
    if complexData
      inoise = rand(dims)-0.5;
    end
  case 'n' % normal distribution with mean 0 and standard dev. 1
    rnoise = randn(dims);
    if complexData
      inoise = randn(dims);
    end
  case 'f' % general 1/f^alpha noise
    alpha = 1;
    rnoise = oneoverfnoise(dims,alpha);
    if complexData
      inoise = oneoverfnoise(dims,alpha);
    end
  otherwise
    error('Wrong noise model. Use ''n'', ''u'', or ''f''.');
end
  
if complexData
  noise = complex(rnoise,inoise);
else
  noise = rnoise;
end

noise = noise/std(noise(:)); % scale to stddev = 1
noise = reshape(noise,dims);

signallevel = max(y(:)) - min(y(:));
signallevel = max([real(signallevel) imag(signallevel)]);

if signallevel>0
  noiselevel = signallevel/SNR;
else
  noiselevel = 1;
end

yn = y + noise*noiselevel;

return
%===============================================================================

% construct frequency domain with 1/f^alpha power characteristics and random phases
function noise = oneoverfnoise(dims,alpha)

nPoints = prod(dims);
freq = 0:floor(nPoints/2);
Amplitudes = freq.^(-alpha/2);
Amplitudes(1) = 1;
Phases = exp(2i*pi*rand(size(Amplitudes)));
FreqDomain = Amplitudes.*Phases;
if mod(nPoints,2)
  d = [FreqDomain conj(FreqDomain(end:-1:2))];
else
  d = [FreqDomain conj(FreqDomain(end-1:-1:2))];
end
noise = real(ifft(d));
isColVec = dims(1)>dims(2);
if isColVec, noise = noise(:); end

return
