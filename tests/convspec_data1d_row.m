function ok = test()

% Test 1D convolution of 1D row data

N = 601;
x = 1:N;
data = zeros(1,N);

% Positions and amplitudes of Gaussian peaks, all with same width
x0 = [205 400];
amp = [1 2.3];
w = 10;

% Build stick spectrum and convolve
data(x0) = amp;
z1 = convspec(data,1,w);

% Build via explicit sum of Gaussians
z2 = zeros(1,N);
for k = 1:numel(x0)
  z2 = z2 + amp(k)*gaussian(x,x0(k),w);
end

ok =  areequal(z1,z2,1e-10,'rel');
