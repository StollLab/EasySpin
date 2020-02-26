function ok = test()

J = 10e3;

S1 = 5/2;
S2 = 2;

Sys.S = [S1 S2];
Sys.ee = J;

[C,E] = spinladder(Sys);

S  = abs(S1-S2):S1+S2;
E0 = J/2*(S.*(S+1)-S1*(S1+1)-S2*(S2+1));

ok = areequal(E,E0,J*1e-10,'abs');
