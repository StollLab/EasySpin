%  autocorrfft  Calculate autocorrelation function using the FFT.
%
%  AutoCorr = autocorrfft(y);
%
%  Input:
%     y              array of data
%     vector         1: treat input data as if the first dimension
%                       represents different components of a vector and use
%                       the dot product (sum over components at the end)
%                    0: treat each component in first dimension separately
%     normalized     1: normalize by the variance (the first data point)
%                    0: no normalization
%     centered       1: subtract the mean before performing FFT
%                    0: no mean subtraction
%
%  Output:
%     autocorr       array

function AutoCorr = autocorrfft(y, vector, normalized, centered)

if nargin==2
  % normalized and centered output by default
  normalized = 1;
  centered = 1;
end

if numel(size(y))>2
  error('The input array must have 1 or 2 dimensions.')
end

if centered
  y = bsxfun(@minus, y , mean(y,2));
end

N = size(y,2);
F = fft(y, 2*N, 2);
r = ifft(F.*conj(F),[],2);
r = real(r(:,1:N));

if vector
  AutoCorr = sum(r,1);
else
  AutoCorr = r;
end

n = N*ones(1, N) - [1:N] + 1;

AutoCorr = AutoCorr./n;

if normalized
  AutoCorr = AutoCorr./AutoCorr(:,1);
end

end