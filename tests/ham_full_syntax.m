function ok = test()

sys = struct('S',1/2,'g',[2 3 4]);
B = rand(1,3)*400;
[H0,Gx,Gy,Gz] = ham(sys);
[H0,G] = ham(sys,B);
H = ham(sys,B);

ok = true;
