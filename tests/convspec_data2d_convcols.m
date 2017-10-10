function [err,data] = test(opt,olddata)

% Test convolution of columns in 2D data
%========================================================
N = 101;
M = 10;
data = zeros(N,M);

% Positions and amplitudes of Gaussian peaks, all with same width
x0 = [25 60];
y0 = [3 7];
amp = [1 2.3];
fwhm = 10;

% Build stick spectrum and convolve
for k = 1:numel(x0)
  data(x0(k),y0(k)) = amp(k);
end
z1 = convspec(data,1,[fwhm 0]);

% Build via explicit sum of Gaussians
z2 = zeros(N,M);
for k = 1:numel(x0)
  c = y0(k);
  z2(:,c) = z2(:,c) + amp(k)*gaussian((1:N).',x0(k),fwhm);
end

err =  ~areequal(z1,z2);
data = [];
