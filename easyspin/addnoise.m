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
%   standard deviation of the noise amplitude distribution. If y is
%   complex-valued, the signal amplitude and the noise are determined
%   separately for the real and imaginary parts, using the larger of
%   the two signal amplitudes to set the noise level for both parts.
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

dims = size(y);

switch noisemodel
  case 'u' % uniform distribution [-0.5,0.5]
    noisefun = @() rand(dims)-0.5;
  case 'n' % normal distribution with mean 0 and standard dev. 1
    noisefun = @() randn(dims);
  case 'f' % general 1/f^alpha noise
    alpha = 1;
    noisefun = @() oneoverfnoise(dims,alpha);
  otherwise
    error('Unknown noise type ''%s''. Use ''n'', ''u'', or ''f''.',noisemodel);
end
  
complexData = ~isreal(y); % check for complex data
if complexData
  noiseRe = noisefun();
  noiseRe = noiseRe/std(noiseRe(:)); % scale to stddev = 1
  noiseIm = noisefun();
  noiseIm = noiseIm/std(noiseIm(:)); % scale to stddev = 1
  noise = complex(noiseRe,noiseIm);
else
  noise = noisefun();
  noise = noise/std(noise(:)); % scale to stddev = 1
end

noise = reshape(noise,dims);

signallevel_re = max(real(y(:))) - min(real(y(:)));
signallevel_im = max(imag(y(:))) - min(imag(y(:)));
signallevel = max(signallevel_re,signallevel_im);

if signallevel>0
  noiselevel = signallevel/SNR;
else
  noiselevel = 1;
end

yn = y + noise*noiselevel;

end
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

end
