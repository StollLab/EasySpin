function ok = test()

% Generate oscillatory signal with offset in imaginary part
t = linspace(0,3,1001);  % µs
f = 2.445;  % MHz
offset = 20i;
V = cos(2*pi*f*t) + offset;

% Phase correct with and without ignoring offset in the imaginary part
[~,ph_noignore] = phasecorr(V,false);
[~,ph_ignore] = phasecorr(V,true);

thresh = 1e-2;
ok(1) = areequal(ph_noignore,-pi/2,thresh,'abs');
ok(2) = areequal(ph_ignore,0,thresh,'abs');

end
