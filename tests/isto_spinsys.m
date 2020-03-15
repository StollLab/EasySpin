function ok = test()

% Test whether isto() handles both spin vectors and spin system structures

Sys.S = [3/2 1];

A = isto(Sys,[3 -2 1]);
B = isto(Sys.S,[3 -2 1]);

ok = areequal(A,B,1e-10,'rel');
