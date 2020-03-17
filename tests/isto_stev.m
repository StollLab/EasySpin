function ok = test()

% Test whether isto() and stev() are consistent

S = 3;
k = 2;
q = 1;

c = -1/2;
T = isto(S,[k,q]);
A = c*(T+T');
B = stev(S,[k,q]);

ok = areequal(A,B,1e-10,'rel');
