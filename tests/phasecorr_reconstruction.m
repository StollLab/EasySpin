function ok = test()

% Vph must equal V.*exp(1i*ph) exactly, for both ignoreOffset settings.

t = linspace(0,3,1001);  % µs
f = 2.445;  % MHz
V = cos(2*pi*f*t).*exp(1i*deg2rad(37)) + 0.3i;

[Vph1,ph1] = phasecorr(V);
ok(1) = areequal(Vph1,V.*exp(1i*ph1),1e-12,'abs');

[Vph2,ph2] = phasecorr(V,true);
ok(2) = areequal(Vph2,V.*exp(1i*ph2),1e-12,'abs');

end
