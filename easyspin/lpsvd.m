% Linear Prediction Singular Value Decomposition
%
%  predictedSpectrum = lpsvd(Spectrum, Time, Method, Order)
% [predictedSpectrum, PredictionParameters] = lpsvd(...)
%
% Performs Linear Prediction SVD using a damped exponential model:
% y = amp * exp(1i*phase) * exp(time * (1i*2*pi*freq - damp) );
%
% Spectrum - The spectrum to be predicted
% Time - The corresponding time vector for the spectrum, necessary to
%        predict the proper phase of the signal.
%
% Method - the method string input determines the LPSVD algorithm used,
%          if not provided it will default to 'ss'
%
% methods:
% 'kt' Based on:
% Kumaresan, R.;Tufts, D.W.; IEEE Trans. Acoust. Speech Signal ASSP-30 833 (1982)
%
% 'ss'  state-space method:
% Kung, S.Y.; Arun, K.S.; Bhaskar Rao, D.V.; J. Opt. Soc. Am. 73, 1799 (1983)
% Barkhuijsen, H.; De Beer, R.; Van Ormondt, D.; J. Mag. Reson. 73, 553 (1987)
%
% 'tls' Hankel Total least squares method:
% Van Huffel, S.;Chen, H.; Decanniere, C.; Van Hecke, P.; J. Mag. Reson.,
% Series A 110, 228 (1994)
%
%
% Order - The number of sinusoides to attempt to fit to the data if no order
%         is provided it will be automatically estimated.
%
% Model Order estimated as per:
% Wax, M.; Kailath, T.; IEEE Trans. Acoust. Speech Signal ASSP-39 387 (1985)
% 'mdl'  minimum description length
% 'aic' Akaike information protocol
% however these methods are known to underestimate the number of components




function [y,parameters] = lpsvd(data,time,method,order)

% Check inputs+

dim = size(data);

if min(dim)>1
  error('Data must be a row or column vector.');
end
if numel(dim)>2
  error('Data must be a row or column vector.');
end

data = data(:);

if length(data)~=length(time)
  error('Time vector must be the same length as data')
end
if nargin<3 || isempty(method)
  method = 'ss';
end

if nargin<4
  order = 'mdl';
else
  m = order;
end

dat = data/max(abs(data));
N = length(data);

% estimate L such that L > M(the number of sinusoidal signals)
L = floor(0.6*N);

% generate the LxM Hankel matrix
A = hankel(conj(dat(2:N-L+1)), conj(dat((N-L+1):N)));

[U,S,V] = svd(A,'econ');
S = diag(S);

% determine the model order
if ischar(order)
  switch order
    case 'aic'
      M = length(S);
      aic = zeros(1,M);
      for k = 0:M-1;
        aic(k+1) = 2*N*( (M-k)*log((sum(S(k+1:M))/(M-k))) - sum(log(S(k+1:M)))) ...
          + 2*k*(2*M-k);
      end
      [dummy, m] = min(aic);
      m = m - 1;
    case 'mdl'
      M = length(S);
      mdl = zeros(1,M);
      for k = 0:M-1;
        mdl(k+1) = N*( (M-k)*log((sum(S(k+1:M))/(M-k))) - sum(log(S(k+1:M))))...
          + k*(2*M-k)*log(N)/2;
      end
      [dummy, m] = min(mdl);
      m = m - 1;
  end
end

% truncate the singular values and eigenvectors to the model order
Um = U(:,1:m);
Sm = S(1:m);
Vm = V(:,1:m);

switch method
  % some optimization can still be done in calculating the signal poles,
  % but for now we're going on the premise that processing power is cheap.
  case 'kt'
    
    h = -conj(dat(1:N-L));
    % solve for the signal poles
    b = Vm*(diag(1./Sm)*(Um'*h));
    s = roots([b(end:-1:1);1]);
    
    % and throw away those that lie inside the unit circle
    s = s(abs(s)<1);
    
  case 'ss'
    
    % remove the top and bottom rows from Um
    Umt = Um(2:end,:);
    Umb = Um(1:end-1,:);
    
    % calculate Zprime
    Zp = (Umb'*Umb)\(Umb'*Umt);
    % diagonalize to yield the signal poles
    s = eig(Zp);
    
  case 'tls'
    % remove the top and bottom rows from Um
    Umt = Um(2:end,:);
    Umb = Um(1:end-1,:);
    
    % obtain the SVD of the augmented matrix
    [dummy,dummy,Vu] = svd([Umb Umt],'econ');
    
    % calculate Zprime
    Zp = -Vu(1:m,m+1:2*m)/Vu(m+1:2*m,m+1:2*m);
    s = eig(Zp);
end
if isempty(s), error('No prediction coefficients calculated, algorithm failed to converge.');end
% pull out the dampings and frequencies
s = -log(s);
dt = time(2) - time(1);
damp = real(s)/dt;
freq = imag(s)/(2*pi*dt);

% reject negative damping
I = damp>0;
damp = damp(I);
freq = freq(I);

% sort the frequencies
[freq,I] = sort(freq);
damp = damp(I);


% generate the signals
y = exp(time(:)*(-damp' + 1i*2*pi*freq'));

% Calculate their phase and amplitude
%a = (y'*y)\(y'*data(:));
% using Cholesky decomposition in hopes of avoiding singular matrices on inversion
R = chol(y'*y);
a = R\(R'\(y'*data(:)));

y = y*a;
amp = abs(a);
phase = imag(log(a./amp));

% attempt to correct and account for negative amplitude
ind = phase<0;
phase(ind) = phase(ind)+pi;
amp(ind) = -amp(ind);

% output
parameters.damping = damp;
parameters.frequency = freq;
parameters.amplitude = amp;
parameters.phase = phase;
parameters.model = @(time) exp(time(:)*(-damp' + 1i*2*pi*freq')) * (amp .*exp(1i*phase));

y = reshape(y,dim);
return










