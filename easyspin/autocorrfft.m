%  autocorrfft  Calculate autocorrelation function using the FFT.
%
%  AutoCorr = autocorrfft(y);
%
%  Input:
%     y              array of data
%
%  Output:
%     autocorr       array

function AutoCorr = autocorrfft(y)

if numel(size(y))>2
  error('The input array must have 1 or 2 dimensions.')
end

  y = bsxfun(@minus, y , mean(y,2));
  N = length(y);
  F = fft(y, 2*N, 2);
  r = ifft(F.*conj(F),[],2);
  r = real(r(:,1:N));
  AutoCorr = sum(r,1).';

n = N*ones(1, N) - [1:N] + 1;

AutoCorr = AutoCorr./n.';
AutoCorr = AutoCorr./AutoCorr(1);

end