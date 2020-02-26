function ok = test()

% Test value

a = hbar;
b = planck/2/pi;
ok = areequal(a,b,1e-12,'rel');
