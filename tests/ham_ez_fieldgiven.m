function ok = test()

% Test two forms of calls

Sys.S = 3/2;
Sys.g = diag([2 2 2])+rand/2;
B = ang2vec(rand*2*pi,rand*pi)*340;

[Gx0,Gy0,Gz0] = ham_ez(Sys);
H0 = Gx0*B(1) + Gy0*B(2) + Gz0*B(3);
H1 = ham_ez(Sys,B);

ok = areequal(H0,H1,1e-10,'abs');
