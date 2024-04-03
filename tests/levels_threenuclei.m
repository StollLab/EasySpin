function ok = test()

Sys.S = 1/2;
Sys.g = [2 3 4];
Sys.Nucs = '14N,14N,63Cu';
Sys.A = [4 7 1; 10 23 -8; 78 78 -300];

phi = 0.4564;  % rad
theta = 1.14343564;  % rad
B = 0;  % mT

E = levels(Sys,[phi theta],B);

ok = true;
