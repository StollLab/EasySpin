% addnoise   Add noise to a signal
%
%   yn = addnoise(y,SNR)
%   yn = addnoise(y,SNR,'f')
%   yn = addnoise(y,SNR,'n')
%   yn = addnoise(y,SNR,'u')
%
%   Adds noise to y to give yn with signal-to-noise ratio SNR.
%   'f': 1/f noise (default)
%   'n': normal (Gaussian) distribution
%   'u': uniform distribution
%
%   The SNR is defined as the ratio of signal amplitude to the
%   standard deviation of the noise amplitude distribution.
%
%   Example:
%     x = linspace(-1,1,1001);
%     y = lorentzian(x,0,0.3,1);
%     yn = addnoise(y,20);
%     plot(x,yn);

function yn = addnoise(y,SNR,model)

if (nargin==0), help(mfilename); return; end

if (nargin<3), model = 'f'; end

if numel(SNR)~=1 || (SNR<=0)
  error('The second input, the signal-to-noise ratio (SNR), must be a positive number.');
end

%rand('state',sum(100*clock)); % obsoleted in R2011a

dims = size(y);
% check for complex data
iscpx = double(~isreal(y));
for cpx=0:iscpx
  if cpx, rnoise = noise;end
  switch model
    case 'u',
      noise = rand(dims)-0.5;
      
    case 'n',
      noise = randn(dims);
      
    case 'f', % general 1/f^alpha noise
      alpha = 1;
      
      % construct frequency domain with 1/f^alpha power characteristics
      % and random phases
      nPoints = numel(y);
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
      
    otherwise
      error('Wrong noise model. Use ''n'', ''u'', or ''f''.');
  end
  if cpx, noise = complex(rnoise,noise);end
end
noise = noise/std(noise(:)); % scale to stddev = 1;
noise = reshape(noise,size(y));

signallevel = (max(y(:))-min(y(:)));
signallevel = max([real(signallevel) imag(signallevel)]);

if (signallevel>0)
  noiselevel = signallevel/SNR;
else
  noiselevel = 1;
end

yn = y + noise*noiselevel;