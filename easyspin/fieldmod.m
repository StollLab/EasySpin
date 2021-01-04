% fieldmod  field modulation 
%
%   yMod = fieldmod(x,y,ModAmpl);
%   yMod = fieldmod(x,y,ModAmpl,Harmonic);
%   fieldmod(...)
%
%   Computes the effect of field modulation
%   on an EPR absorption spectrum.
%
%   Input:
%   - x: magnetic field axis vector [mT]
%   - y: absorption spectrum
%   - ModAmpl: peak-to-peak modulation amplitude [mT]
%   - Harmonic: harmonic (0, 1, 2, ...); default is 1
%
%   Output:
%   - yMod: pseudo-modulated spectrum
%
%   If no output variable is given, fieldmod plots the
%   original and the modulated spectrum.
%
%   Example:
%
%     x = linspace(300,400,1001);
%     y = lorentzian(x,342,4);
%     fieldmod(x,y,20);

% References
% --------------------------------------------------
% Berger, Günthart, Z.Angew.Math.Phys. 13, 310 (1962)
% Wilson, J.Appl.Phys. 34, 3276 (1963)
% Haworth, Richards, Prog.Nmr.Spectrosc. 1, 1 (1966)
% Hyde et al., Appl.Magn.Reson. 1, 483-496 (1990)
% Hyde et al., J.Magn.Reson. 96, 1-13 (1992)
% Kaelin, Schweiger, J.Magn.Reson. 160, 166-180 (2003)
% Nielsen, Robinson, Conc. Magn. Reson. A 23, 38-48 (2004)

function varargout = fieldmod(x,y,ModAmpl,Harmonic)

if nargin==0, help(mfilename); return; end

Display = (nargout==0);

if nargin<3 || nargin>4, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>1, error('Too many output arguments.'); end

% Supplement arguments and check range.
if nargin<4, Harmonic = 1; end
if numel(Harmonic)~=1 || Harmonic<0 || ~isreal(Harmonic) || mod(Harmonic,1)
  error('Harmonic must be a positive integer (1, 2, 3, etc)!');
end

% Check ModAmpl
if ModAmpl<=0
  error('Modulation amplitude (3rd argument) must be positive.');
end

% Get length of vectors.
n = length(x);
if length(y)~=n, error('x and y must have the same length!'); end

sizey = size(y);
if all(sizey~=1)
  error('y must be a row or column vector.');
end

isRowVector = (sizey(1)==1);
y = y(:);

% Compute relative base-to-peak amplitude.
dx = x(2) - x(1);
Ampl = ModAmpl/2/dx;

% FFT-based convolution
%------------------------------------------------------------
% Compute FFT of input signal, zero negative part.
NN = 2*n+1; % to avoid fold-around during convolution
ffty = fft(y,NN);
ffty(ceil(NN/2)+1:end) = 0;

% Convolution with IFT of Bessel function.
S = (0:NN-1).'/NN;
yMod = ifft(ffty.*besselj(Harmonic,2*pi*Ampl*S));
yMod = yMod(1:n); % pick out the correct subarray

% Adjust phase.
yMod = (1i)^Harmonic * yMod;

if isRowVector
  yMod = yMod.';
end

yModInPhase = real(yMod);
%yModOutOfPhase = imag(yMod);

if Display
  clf
  subplot(2,1,1);
  plot(x,y);
  xlim([min(x) max(x)]);
  title('Original spectrum');
  subplot(2,1,2);
  plot(x,yModInPhase);
  xlim([min(x) max(x)]);
  xlabel('magnetic field (mT)');
  title(sprintf('Modulated spectrum, harmonic %d, modulation amplitude %g mT',Harmonic,ModAmpl));
end

if nargout==1
  varargout = {yModInPhase};
end

return
