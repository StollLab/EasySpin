function [err,data] = test(opt,olddata)

% Test 1D convolution of 1D column data
%========================================================
N = 601;
x = (1:N).';
data = zeros(N,1);

% Positions and amplitudes of Gaussian peaks, all with same width
x0 = [205 456];
amp = [1 2.3];
w = 15;

% Build stick spectrum and convolve
data(x0) = amp;
z1 = convspec(data,1,w);

% Build via explicit sum of Gaussians
z2 = zeros(N,1);
for k = 1:numel(x0)
  z2 = z2 + amp(k)*gaussian(x,x0(k),w);
end

err =  ~all(size(z1)==size(z2)) || ~areequal(z1,z2,1e-10,'rel');
data = [];
