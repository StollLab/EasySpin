function ok = test()

% Test Bkq terms

rng(454);

S = 7/2;
B43 = rand;
B40 = rand;

Sys.S = S;
Sys.B4 = [0 B43 0 0 0 0 0 0 0];
H1 = zfield(Sys);
H0 = B43*stev(S,[4,3]);

ok(1) = areequal(H1,H0,1e-10,'abs');

Sys.B4 = [0 0 0 0 B40 0 0 0 0];
H1 = zfield(Sys);
H0 = B40*stev(S,[4,0]);

ok(2) = areequal(H1,H0,1e-10,'abs');

Sys.B4 = B40;
H1 = zfield(Sys);
H0 = B40*stev(S,[4,0]);
ok(3) = areequal(H1,H0,1e-10,'abs');


