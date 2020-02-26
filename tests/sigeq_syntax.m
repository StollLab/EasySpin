function ok = test()

[Sx,Sy,Sz] = sop([1/2 1],'xe','ye','ze');
H = Sx + 5*Sy - 2*Sz;

sigeq(H,10);
sigeq(H,10,'pol');
sigma = sigeq(H,10);
sigma = sigeq(H,10,'pol');

ok = true;
