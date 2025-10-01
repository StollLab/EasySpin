function ok = test()

% Test two-dimensional polynomial baseline correction

x = linspace(0,1,7).';
y = linspace(0,1,6);
data = 1 + 0.12*x - 0.9*y + x.^2 + 0.445*y.^2 - 0.2*x.*y;

dim = [];
n = [2 2];
datacorr = basecorr(data,dim,n);

ok = max(abs(datacorr))<1e-10;