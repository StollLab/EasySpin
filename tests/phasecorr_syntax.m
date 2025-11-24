function ok = test()

t = linspace(0,3,1001);  % µs
f = 2.445;  % MHz
offset = 0.1i;
V = cos(2*pi*f*t) + offset;

Vph = phasecorr(V);
[Vph,ph] = phasecorr(V);
Vph = phasecorr(V,true);
Vph = phasecorr(V,false);

ok = true;

end
