function ok = test()

J = 50*30e3;

Sys.S = [3/2 3/2];

Sys.D = [10 1; 20 6];
Sys.ee = rand(3)*30e3 + J;
F = spinladder(Sys);

Sys.D = rand(2,3)*10;
Sys.ee = rand(3)*30e3 + J;
F = spinladder(Sys);

Sys.D = [10 20];
Sys.ee = rand(3)*30e3 + J;
F = spinladder(Sys);

ok = true;
