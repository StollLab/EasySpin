% rapidscan2spc  Conversion of sinusoidal rapid-scan signal to EPR spectrum
%
%   [B,spc] = rapidscan2spc(rsSignal,rsAmp,rsFreq,g)
%
% Converts a sinusoidal rapid-scan time-domain EPR signal to the corresponding
% field-swept EPR spectrum.
%
% Inputs:
%   rsSignal time-domain rapid-scan signal (My+1i*Mx) over one period of
%            the field modulation (vector); My is the absorption, and Mx
%            is the dispersion
%   rsAmp    peak-to-peak modulation amplitude (in mT)
%   rsFreq   modulation frequency (in kHz)
%   g        (optional) g value, needed for field-domain rapid scan.
%            If omitted, g = 2.0023193 is assumed.
%
% Outputs:
%   dB       field offset axis (in mT)
%   spc      spectrum (both absorption and dispersion)
%
% For this deconvolution to work, the rapid-scan signals of the up- and
% the down-sweeps should not overlap, and the data should not be saturated.
%
% Additionally, the input signal needs to be correctly quadrature phased, time
% shifted to align with a cosine field modulation, and background corrected.

function [dB,spc] = rapidscan2spc(M,rsAmp,rsFreq,g)

%{
This function implements the deconvolution method from
  Tseitlin, Rinard, Quine, Eaton, Eaton
  Deconvolution of sinusoidal rapid EPR scans
  J. Magn. Reson. 208, 279-283, 2011
  https://doi.org/10.1016/j.jmr.2010.11.015
%}

if nargin==0
  help(mfilename);
  return
end

if nargin<3
  error('At least three inputs are expected: M, rsAmp, and rsFreq');
end

if nargin<4
  g = gfree;
end

if ~isnumeric(M) || ~isvector(M)
  error('First input (rapid-scan signal) must be a vector.');
end

if isreal(M)
  error('First input (rapid-scan signal) must contain both absorption and dispersion.');
end

M = M(:);
N = numel(M);
if mod(N,2)
  error('First input (rapid-scan signal) must contain an even number of points.');
end

if numel(rsAmp)~=1 || ~isreal(rsAmp) || rsAmp<=0
  error('Second input (modulation amplitude) must be a positive number.')
end

if numel(rsFreq)~=1 || ~isreal(rsFreq) || rsFreq<=0
  error('Third input (modulation frequency) must be a positive number.')
end

rsFreq_kHz = rsFreq;
rsFreq = rsFreq_kHz*1e3; % kHz -> Hz

% Construct time vector over one period
Tperiod = 1/rsFreq; % s
t = (0:N-1).'/N*Tperiod; % s

% Convert modulation amplitude to rad/s
rsAmp_mT = rsAmp;
rsAmp = rsAmp_mT/1e3; % mT -> T
A = rsAmp/2;
gamma = g*bmagn/hbar;
A = gamma*A; % -> rad/s

% Phase transformation
om = 2*pi*rsFreq;
ph = A*sin(om*t)/om; % = integral of A*cos(omMod*t)
phf = exp(-1i*ph);
Mph = M.*phf;

% Fourier deconvolution of the two half-periods and add up
fourierdeconv = @(idx) fftshift(fft(Mph(idx))./fft(phf(idx)));
spc = fourierdeconv(1:N/2) + fourierdeconv(N/2+1:N);

% Generate field axis and truncate to mod. amplitude range
nu = fdaxis(t(1:N/2))/1e6; % Hz -> MHz
dB = -mhz2mt(nu,g);
idx = abs(dB)<=rsAmp_mT/2;
dB = dB(idx);
spc = spc(idx);

% Plotting
%-------------------------------------------------------------------------------
if nargout==0
  
  subplot(2,1,1);
  t_us = t*1e6;
  plot(t_us,real(M),t_us,imag(M));
  grid on
  xlabel('time (us)');
  title('time-domain signal');
  legend('absorption','dispersion');
  legend boxoff
  
  subplot(2,1,2);
  plot(dB,real(spc),dB,imag(spc));
  grid on
  xlabel('field offset (mT)');
  title('spectrum');
  legend('absorption','dispersion');
  legend boxoff

end
