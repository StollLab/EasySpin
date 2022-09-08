function ok = test()

Sys = struct('S',[1/2 1/2 1/2],'g',[2 2 2; 2.1 2.1 2.1; 2.2 2.2 2.2]);
Sys.ee = [1 1 1; 2 2 2; 3 3 3];
H = ham_ee(Sys);
H = ham_ee(Sys,[1 2]);
H = ham_ee(Sys,[1 3]);
H = ham_ee(Sys,[2 3]);

ok = true;
