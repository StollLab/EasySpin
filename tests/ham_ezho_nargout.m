function ok = test()

% Assert that ham_ezho returns as many tensor orders as outputs are requested.

rng(5);

Sys.S = 1;
Sys.Ham112 = rand(1,5);
Sys.Ham312 = rand(1,5);

[G0,G1,G2] = ham_ezho(Sys);
cG = ham_ezho(Sys);

threshold = 1e-10;
ok(1) = numel(cG)==4;
ok(2) = areequal(G0,cG{1},threshold,'abs');
ok(3) = areequal(G1{2},cG{2}{2},threshold,'abs');
ok(4) = areequal(G2{1,3},cG{3}{1,3},threshold,'abs');
