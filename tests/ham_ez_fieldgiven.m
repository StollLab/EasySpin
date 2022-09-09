function ok = test()

% Test two forms of calls

Sys.S = 3/2;
Sys.g = diag([2 2 2])+rand/2;
B = ang2vec(rand*2*pi,rand*pi)*340;

[mux0,muy0,muz0] = ham_ez(Sys);
H0 = -(mux0*B(1) + muy0*B(2) + muz0*B(3));
H1 = ham_ez(Sys,B);

ok = areequal(H0,H1,1e-10,'abs');
