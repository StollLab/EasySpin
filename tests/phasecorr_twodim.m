function ok = test()

t = linspace(0,3,1001);  % µs
f = 2.445;  % MHz
V = cos(2*pi*f*t).*exp(1i*deg2rad(12));

% Generate 2D array
V = V + 0.1*randn(30,length(t));

Vph = phasecorr(V);

ok = ismatrix(Vph);

end
