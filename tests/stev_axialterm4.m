function ok = test()

% Compare two expressions for 4th-order axial term

S = 5/2;

[E,Sx,Sy,Sz] = sop(S,'e','x','y','z');

s = S*(S+1);
P1 = 35*Sz^4 - (30*s - 25)*Sz^2 - E*(6*s - 3*s^2);

P2 = stev(S,[4,0]);

ok = areequal(P1,P2,1e-10,'rel');
