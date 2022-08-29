function ok = test()

sys = struct('S',1/2,'g',[2 3 4]);
B = rand(1,3)*400;
[F,Gx,Gy,Gz] = ham(sys);
[F,G] = ham(sys,B);
H = ham(sys,B);

ok = true;
