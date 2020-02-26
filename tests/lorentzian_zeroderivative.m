function ok = test()

% Zero derivative at center

a = lorentzian(0,0,1.234,1);
ok = abs(a)<1e-15;
