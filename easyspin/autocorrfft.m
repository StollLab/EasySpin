%  autocorrfft  Calculate autocorrelation function using the FFT.
%
%  AutoCorr = autocorrfft(y);
%
%  Input:
%     y              array of data
%
%  Output:
%     autocorr       array

function AutoCorr = autocorrfft(y);

if numel(size(y)) ~= 2
  error('The input array must be 2-dimensional.')
end

if size(y, 1) > 1 & size(y, 2) > 1
  error('The input array must a vector.')
end

y = bsxfun(@minus, y , mean(y));

N = length(y);

r = ifft(fft(y, 2*N).*conj(fft(y, 2*N)));
AutoCorr = real(r(1:N));

n = N*ones(1, N) - [1:N] + 1;

AutoCorr = AutoCorr./n;
AutoCorr = AutoCorr./AutoCorr(1);

end