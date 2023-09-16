function ok = test()

% Test two-dimensional polynomial baseline correction

x = linspace(0,5,7).';
y = linspace(0,10,6);
data = 1 + 0.12*x - 0.9*y + x.^2 + y.^2 - 0.2*x.*y;

corrdata = basecorr(data,[],[2 2]);

ok = max(abs(corrdata))<1e-10;
