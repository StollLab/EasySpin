% fieldmod  field modulation 
%
%   spcMod = fieldmod(B,spc,ppModAmp);
%   spcMod = fieldmod(B,spc,ppModAmp,Harmonic);
%   spcMod = fieldmod(B,spc,ppModAmp,Harmonic,Method);
%   fieldmod(...)
%
%   Computes the effect of field modulation on an EPR absorption spectrum.
%
%   Input:
%   - B: magnetic field axis vector, mT
%   - spc: absorption spectrum
%   - ppModAmp: peak-to-peak modulation amplitude, mT
%   - Harmonic: harmonic (0, 1, 2, ...); default is 1
%   - Method: calculation method, 'fft' or 'conv'; default is 'fft'
%
%   Output:
%   - spcMod: pseudo-modulated spectrum
%
%   If no output variable is given, fieldmod plots the original and
%   the modulated spectrum.
%
%   Example:
%
%     B = linspace(300,400,1001);  % mT
%     spc = lorentzian(B,342,4);
%     fieldmod(B,spc,20);

% The FFT method is generally more accurate, its only downside is that
% it causes Gibbs wiggles in certain cases.

% References
% --------------------------------------------------
% Berger, Günthart, Z.Angew.Math.Phys. 13, 310 (1962)
% Wilson, J.Appl.Phys. 34, 3276 (1963)
% Haworth, Richards, Prog.Nmr.Spectrosc. 1, 1 (1966)
% Hyde et al., Appl.Magn.Reson. 1, 483-496 (1990)
% Hyde et al., J.Magn.Reson. 96, 1-13 (1992)
% Kaelin, Schweiger, J.Magn.Reson. 160, 166-180 (2003)
% Nielsen, Robinson, Conc. Magn. Reson. A 23, 38-48 (2004)

function varargout = fieldmod(B,spc,ppModAmp,Harmonic,Method)

if nargin==0, help(mfilename); return; end

% Check input and output arguments
if nargin<3 || nargin>5, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>1, error('Too many output arguments.'); end

plotResult = nargout==0;

% Supplement arguments and check range
if nargin<4, Harmonic = 1; end
if numel(Harmonic)~=1 || Harmonic<0 || ~isreal(Harmonic) || mod(Harmonic,1)
  error('Harmonic must be a positive integer (1, 2, 3, etc)!');
end
if nargin<5, Method = 'fft'; end

% Check ppModAmp
if ppModAmp<=0
  error('Peak-to-peak modulation amplitude (3rd argument) must be positive.');
end

% Check Method
if Method~="fft" && Method~="conv"
  error(' Method (5th input) must be either ''fft'' or ''conv''.');
end

% Get length of vectors
n = length(B);
if length(spc)~=n, error('x and y must have the same length!'); end

if ~isvector(spc)
  error('spc (2nd input) must be a row or column vector.');
end

isRowVector = isrow(spc);
spc = spc(:);

dB = mean(diff(B));  % field increment

switch Method
  case 'fft'
    % Convolution with Bessel function via FFT
    %---------------------------------------------------------------------------
    % Compute FFT of input signal, zero negative part
    NN = 2*n+1; % to avoid fold-around during convolution
    spc_fft = fft(spc,NN);
    spc_fft(ceil(NN/2)+1:end) = 0;

    % Multiply with Bessel function and inverse FT
    S = (0:NN-1).'/NN;
    bpAmp = ppModAmp/2/dB;  % base-to-peak amplitude relative to increment
    spc_fft = spc_fft.*besselj(Harmonic,2*pi*bpAmp*S);
    yMod = ifft(spc_fft);
    yMod = yMod(1:n);  % pick out the positive subarray
    yMod = (1i)^Harmonic * yMod;  % adjust phase

  case 'conv'
    % Direct convolution with modulation kernel in field domain
    %---------------------------------------------------------------------------
    % This method is less accurate than the Fourier transform, therefore we
    % interpolate here to improve accuracy. Its only advantage over the
    % Fourier transform method is that it avoids Gibbs wiggles when the
    % spectrum is much narrower than the modulation amplitude, or if the
    % spectrum is non-zero at the edges of the field/frequency range.
    Ni = 10; % interpolation factor
    spc = interp1(spc,1:1/Ni:length(spc),'makima');  % cubic interpolation
    dx = dB/ppModAmp/Ni;
    x = -1+dx:2*dx:1-dx;
    modkernel = (-1)^Harmonic*chebyT(Harmonic,x)./sqrt(1-x.^2);
    yMod = conv(spc,modkernel,'same');
    yMod = yMod(1:Ni:end);  % downsample again
    yMod = yMod*dx/pi;

end

if isRowVector
  yMod = yMod.';
end

yModInPhase = real(yMod);
%yModOutOfPhase = imag(yMod);

% Plotting
%-------------------------------------------------------------------------------
if plotResult
  clf
  subplot(2,1,1);
  plot(B,spc);
  xlim([min(B) max(B)]);
  title('Original spectrum');
  subplot(2,1,2);
  plot(B,yModInPhase);
  xlim([min(B) max(B)]);
  xlabel('magnetic field (mT)');
  title(sprintf('Modulated spectrum, harmonic %d, modulation amplitude %g mT',Harmonic,ppModAmp));
end

if nargout==1
  varargout = {yModInPhase};
end

end

% chebyT Chebyshev polynomial of the first kind
%
%   y = chebyT(n,x)
%
%   Evaluates the Chebyshev polynomial of the first kind T_n(x)
%
%   Input:
%     n   order of Chebyshev polynomial (0, 1, 2, ...)
%     x   array of values to evaluate T_n(x) for
%   Output:
%     y   array of values of T_n(x), same size as x
function y = chebyT(n,x)

if n<0 || mod(n,1)~=0
  error('Chebyshev polynomial order n must be 0, 1, 2, ...');
end

y0 = ones(size(x));
if n==0, y = y0; return; end

y1 = x;
if n==1, y = y1; return; end

% Iterative evaluation using recurrence relation
for k = 2:n
  y = 2*x.*y1 - y0;
  y0 = y1;
  y1 = y;
end

end
