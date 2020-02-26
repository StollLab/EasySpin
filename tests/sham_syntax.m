function ok = test()

sys = struct('S',1/2,'g',[2 3 4]);
B = rand(1,3)*400;
[F,Gx,Gy,Gz] = sham(sys);
[F,G] = sham(sys,B);
H = sham(sys,B);

ok = true;
