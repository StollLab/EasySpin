function ok = test()

Sys.Nucs = '1H,14N';
Sys.A = [12 23 34; 3 5 7];
Sys.Q = [0 0 0; -1 -1 2];

B0 = ang2vec(pi/3,pi/6)*300;  % mT

Hnz1 = ham_nz(Sys,B0,1);
Hnz2 = ham_nz(Sys,B0,2);
Hnz = ham_nz(Sys,B0);

ok = areequal(Hnz,Hnz1+Hnz2,1e-10,'abs');
