function ok = test()

% Test two-dimensional polynomial baseline correction

x = linspace(0,1,201).';
y = linspace(0,1,101);
data0 = gaussian(x,0.4,0.1)*gaussian(y,0.6,0.1);
baseline = 1 + 0.12*x - 0.9*y + x.^2 + 0.445*y.^2 - 0.2*x.*y;

data = data0 + baseline;

dim = [];
n = [2 2];
baselineCheck = basecorr(baseline,dim,n);

regionx = x<0.1 | x>0.7;
regiony = y<0.2 | y>0.9;
region = regionx | regiony;

[datac, baselinec] = basecorr(data, dim, n, region);

ok(1) = max(abs(baselineCheck), [], 'all')<1e-10;
threshold = 1e-10;
ok(2) = areequal(datac,data0,threshold,'abs');
ok(3) = areequal(baseline,baselinec,threshold,'abs');
