function ok = test()

% Compare two expressions for 4th-order cubic term
S = 5/2;

[E,Sx,Sy,Sz] = sop(S,'e','x','y','z');

P1 = Sx^4 + Sy^4 + Sz^4 - E*1/5*S*(S+1)*(3*S^2+3*S-1);
P1 = P1/6;

P2 = stev(S,[4,0]) + 5*stev(S,[4,4]);
P2 = P2/120;

ok = areequal(P1,P2,1e-10,'rel');

