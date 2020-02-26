function ok = test()

Sys = struct('S',1,'g',[3 4 5]);
v = spinvec(Sys);
ok = (v==1);
