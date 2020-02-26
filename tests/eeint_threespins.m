function ok = test()

Sys = struct('S',[1/2 1/2 1/2],'g',[2 2 2; 2.1 2.1 2.1; 2.2 2.2 2.2]);
Sys.ee = [1 1 1; 2 2 2; 3 3 3];

a = [1 2 3 2 3 4 5 6 7 6 7 8];
b = [1 2 2 3 3 4 5 6 6 7 7 8];
c = [1 -1 2 2 -1 1 1 -1 2 2 -1 1]*3/4;
H0_23 = full(sparse(a,b,c));

a = [1 2 5 3 4 7 2 5 6 4 7 8];
b = [1 2 2 3 4 4 5 5 6 7 7 8];
c = [1 -1 2 1 -1 2 2 -1 1 2 -1 1]*1/2;
H0_13 = full(sparse(a,b,c));

a = [1 2 3 5 4 6 3 5 4 6 7 8];
b = [1 2 3 3 4 4 5 5 6 6 7 8];
c = [1 1 -1 2 -1 2 2 -1 2 -1 1 1]/4;
H0_12 = full(sparse(a,b,c));

H1_23 = eeint(Sys,[2 3]);
H1_13 = eeint(Sys,[1 3]);
H1_12 = eeint(Sys,[1 2]);

thr = 1e-10;
ok = areequal(H0_12,H1_12,thr,'abs') && ...
     areequal(H0_13,H1_13,thr,'abs') && ...
     areequal(H0_23,H1_23,thr,'abs');
