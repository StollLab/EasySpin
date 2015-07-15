% hilberttrans  Hilbert transform
%
%   yh = hilberttrans(y)
%
%   Computes the Hilbert transform of input vector y. It can
%   be used to compute a dispersion signal from an absorption
%   signal.
%
%   Input:
%     - y:    input vector (spectrum)
%   Output:
%     - yh:   Hilbert transform
%               real part: identical to y
%               imaginary part: actual Hilbert transform of y
%
%   Example:
%      x = linspace(-1,1,1000);
%      yabs = gaussian(x,0,0.1,0);
%      ydisp = hilberttrans(yabs);
%      ydisp = imag(ydisp);
%      plot(x,yabs,x,ydisp)

function y_h = hilberttrans(y)

if nargin==0
  help(mfilename);
  return;
end

if ~isreal(y)
  error('Input must be real.');
end

dim = size(y);

if min(dim)>1
  error('Input must be a row or column vector.');
end
if numel(dim)>2
  error('Input must be a row or column vector.');
end


% Build Hilbert kernel
n = max(dim);
if mod(n,2)==0
  idx = n/2;
  h(1) = 1;
  h(idx+1) = 1;
  h(2:idx) = 2;
  h(idx+2:n) = 0;
else
  idx = ceil(n/2);
  h(1) = 1;
  h(2:idx) = 2;
  h(idx+1:n) = 0;
end

% Apply via Fourier Transform
y_h = ifft(fft(y(:)).*h(:));

% shift back to original dimensions
y_h = reshape(y_h,dim);

return

