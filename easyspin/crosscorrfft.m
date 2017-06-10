%  crosscorrfft  Calculate cross-correlation function using the FFT.
%
%  CrossCorr = crosscorrfft(y);
%
%  This function assumes that the second dimension is the "long" axis.
%
%  Input:
%     x              array of data
%     y              array of data
%     normalized     1: normalize by the variance (the first data point)
%                    0: no normalization
%     vector         1: treat input data as if the first dimension
%                       represents different components of a vector and use
%                       the dot product (sum over components at the end)
%                    0: treat each component in first dimension separately
%     centered       1: subtract the mean before performing FFT
%                    0: no mean subtraction
%
%  Output:
%     CrossCorr       array

function CrossCorr = autocorrfft(x, y, vector, normalized, centered)

if nargin==3
  % normalized and centered output by default
  normalized = 1;
  centered = 1;
end

if ndims(x)~=2||ndims(y)~=2
  error('The input arrays must have 2 dimensions.')
end

if isequal(size(x),size(y))
  error('Input arrays must be the same size.')
end

if centered
  x = bsxfun(@minus, x , mean(x,2));
  y = bsxfun(@minus, y , mean(y,2));
end

N = length(y);
Fx = fft(x, 2*N, 2);
Fy = fft(y, 2*N, 2);
r = ifft(Fx.*conj(Fy),[],2);
r = real(r(:,1:N));

if vector
  CrossCorr = sum(r,1);
else
  CrossCorr = r;
end

n = N*ones(1, N) - [1:N] + 1;

CrossCorr = CrossCorr./n;

if normalized
  CrossCorr = CrossCorr./CrossCorr(:,1);
end

end