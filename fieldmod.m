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
%   - Harmonic: harmonic, greater than 0, default is 1
%
%   Output;
%   - yMod: pseudo-modulated spectrum
%
%   If no output variable is give, fieldmod plots the
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
% Haworth, Richards, Prog.Nmr.Spec. 1, 1 (1966)
% Hyde et al., Appl.Magn.Reson. 1, 483-496 (1990)
% Hyde et al., J.Magn.Reson. 96, 1-13 (1992)

function varargout = fieldmod(x,y,ModAmpl,Harmonic)

if (nargin==0), help(mfilename); return; end

Display = (nargout==0);

if (nargin<3) || (nargin>4), error('Wrong number of input arguments!'); end
if (nargout<0), error('Not enough output arguments.'); end
if (nargout>1), error('Too many output arguments.'); end

% Supplement arguments and check range.
if (nargin<4), Harmonic = 1; end
if numel(Harmonic)~=1 || (Harmonic<1) || ~isreal(Harmonic) || mod(Harmonic,1)
  error('Harmonic must be a positive integer (1, 2, 3, etc)!');
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

% Compute FFT of input signal, zero negative part.
NN = 2*n+1; % to avoid fold-around during convolution
ffty = fft(y,NN);
ffty(ceil(NN/2)+1:end) = 0;

% Compute relative base-to-peak amplitude.
Ampl = ModAmpl/2/(x(2)-x(1));

% Convolution with IFT of Bessel function.
S = (0:NN-1).'/NN;
yMod = ifft(ffty.*besselj(Harmonic,2*pi*Ampl*S));

% Adjust phase and pick out the correct part.
yMod = (1i)^Harmonic * yMod(1:n);
yModInPhase = real(yMod); % in-phase component
%yModOutPhase = imag(yMod); % out-of-phase component

if isRowVector
  yModInPhase = yModInPhase.';
end

if (Display)
  subplot(3,1,1);
  plot(x,y);
  title('Original spectrum');
  subplot(3,1,[2 3]);
  plot(x,yModInPhase);
  xlabel('magnetic field [mT]');
  title(sprintf('Modulated spectrum, harmonic %d, modulation amplitude %g mT',Harmonic,ModAmpl));
end

if (nargout==1), varargout = {yModInPhase}; end

return
