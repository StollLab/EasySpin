function ok = test()

% Test value

a = hartree;
b = 4.3597447222071e-18;
ok = areequal(a,b,1e-13,'rel');
