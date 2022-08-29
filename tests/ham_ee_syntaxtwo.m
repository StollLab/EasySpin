function ok = test()

Sys = struct('S',[1/2 1/2],'g',[2 2 2;2.1 2.1 2.1],'ee',[10 10 10]);

ham_ee(Sys);
ham_ee(Sys,[1 2]);
H = ham_ee(Sys);
H = ham_ee(Sys,[1 2]);

ok = true;
