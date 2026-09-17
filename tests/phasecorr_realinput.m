function ok = test()

% A purely real signal (no offset) needs no phase correction.

t = linspace(0,3,1001);  % µs
f = 2.445;  % MHz
V = cos(2*pi*f*t);

[Vph,ph] = phasecorr(V);
ok(1) = areequal(ph,0,1e-8,'abs');
ok(2) = areequal(Vph,V,1e-8,'abs');

end
