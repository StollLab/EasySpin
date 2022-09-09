function ok = test()

rng(3454)
Sys = struct('S',1/2,'g',[2 3 4]);
B = rand(1,3)*400;
[H0,mux,muy,muz] = ham(Sys);
[H0,muzL] = ham(Sys,B);
H = ham(Sys,B);

ok = true;
