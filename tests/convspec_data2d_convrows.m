function [err,data] = test(opt,olddata)

% Test convolution of rows in 2D data
%========================================================
N = 10;
M = 101;
data = zeros(N,M);

% Positions and amplitudes of Gaussian peaks, all with same width
x0 = [25 60];
y0 = [3 7];
amp = [1 2.3];
fwhm = 10;

% Build stick spectrum and convolve
for k = 1:numel(x0)
  data(y0(k),x0(k)) = amp(k);
end
z1 = convspec(data,1,[0 fwhm]);

% Build via explicit sum of Gaussians
z2 = zeros(N,M);
for k = 1:numel(x0)
  r = y0(k);
  z2(r,:) = z2(r,:) + amp(k)*gaussian(1:M,x0(k),fwhm);
end

err =  ~areequal(z1,z2);
data = [];
