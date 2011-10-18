% addnoise   Add noise to a signal
%
%   yn = addnoise(y,SNR,'u')
%   yn = addnoise(y,SNR,'n')
%
%   Adds noise to y to give yn with signal-to-noise ratio SNR.
%   'u': uniform distribution between -0.5/SNR and 0.5/SNR
%   'n': normal distribution with standard deviation 1/SNR
%
%   Example:
%     x = linspace(-1,1,1001);
%     y = gaussian(x,0,0.3);
%     yn = addnoise(y,10);
%     plot(x,yn);

function yn = addnoise(y,snr,model)

if (nargin==0), help(mfilename); return; end

if (nargin<3), model = 'u'; end

signallevel = max(y(:))-min(y(:));
noiselevel = signallevel/snr;

%rand('state',sum(100*clock)); % obsoleted in R2011a

switch model
  case 'u', noise = rand(size(y))-0.5;
  case 'n', noise = randn(size(y))*1/2;
  otherwise
    error('Wrong noise model. Use either ''n'' or ''u''.');
end

yn = y + noise*noiselevel;

return
