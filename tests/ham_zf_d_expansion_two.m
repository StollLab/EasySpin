function ok = test()

% Test whether the different input format for D work

% one electron
Sys.S = [3/2 1];
Sys.ee = 1;

D = 100*rand(1,2);


Sys.D = D;
H1 = ham_zf(Sys);
Sys.D = [D(1) 0; D(2), 0];
H2 = ham_zf(Sys);
Sys.D = D(:)*[-1,-1,2]/3;
H3 = ham_zf(Sys);

Sys = rmfield(Sys,'D');
Sys.D_ = [D(1) 0; D(2), 0];
H4 = ham_zf(Sys);

% test
threshold = 1e-7;
ok(1) = areequal(H1,H2,threshold,'abs');
ok(2) = areequal(H1,H3,threshold,'abs');
ok(3) = areequal(H1,H4,threshold,'abs');
