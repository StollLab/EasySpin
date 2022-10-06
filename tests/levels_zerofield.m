function ok = test()

% Energy levels at zero field

Sys.S = 5/2;
Sys.g = 2;
Sys.D = 3e3*[1 0.2];

E = levels(Sys);

E0 = eig(ham(Sys,[0;0;0]));

ok = areequal(E,E0,1e-10,'abs');
