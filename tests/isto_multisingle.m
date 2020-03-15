function ok = test()

S = [2 3/2 1];
k = 3;
q = 1;

kqi = [k q 2];
kq = [0 0; k q; 0 0];

Op = isto(S,kqi,'sparse');
Op0 = isto(S,kq,'sparse');

ok = areequal(Op,Op0,1e-10,'abs');
