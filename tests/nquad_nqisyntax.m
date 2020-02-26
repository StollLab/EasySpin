function ok = test()

% Nuclear quadrupole Hamiltonian matrix

eeqQ = 17;
eta = 0.1;

Sys.Nucs = '2H';
I = nucspin(Sys.Nucs);
Sys.A = 1;

Sys.Q = eeqQ;
H = nquad(Sys);
Ha1 = H(1:3,1:3);
Hb1 = diag([1 -2 1])*eeqQ/(4*I*(2*I-1));

Sys.Q = [eeqQ eta];
H = nquad(Sys);
Ha2 = H(1:3,1:3);
Hb2 = [1 0 eta ; 0 -2 0; eta 0 1]*eeqQ/(4*I*(2*I-1));

ok = areequal(Ha1,Hb1,1e-8,'abs') && areequal(Ha2,Hb2,1e-8,'abs');
