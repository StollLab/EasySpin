function [err,data] = test(opt,olddata)

% Test 2D convolution of 2D data
%========================================================
N = 20;
M = 81;
data = zeros(N,M);

% Positions and amplitudes of Gaussian peaks, all with same width
x0 = [10 14];
y0 = [25 60];
amp = [1 2.3];
fwhmx = 8;
fwhmy = 7;

% Build stick data array
for k = 1:numel(x0)
  data(x0(k),y0(k)) = amp(k);
end

% 2D convolution
z1 = convspec(data,1,[fwhmx fwhmy]);

% Convolution along rows, followed by convolution along columns
z2 = convspec(data,1,[0 fwhmy]);
z2 = convspec(z2,1,[fwhmx 0]);

% Convolution along colums, followed by convolution along rows
z3 = convspec(data,1,[fwhmx 0]);
z3 = convspec(z3,1,[0 fwhmy]);

err =  ~areequal(z1,z2) || ~areequal(z1,z3);
data = [];
