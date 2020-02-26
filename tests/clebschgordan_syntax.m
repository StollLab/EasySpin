function ok = test()

a1 = clebschgordan(4,3,2,1,-2,-1);
a2 = clebschgordan([4 3 2],[1 -2 -1]);
a3 = clebschgordan([4 1],[3 -2],[2 -1]);
a4 = -sqrt(7)/6;

ok = areequal([a1 a2 a3],[a4 a4 a4],1e-10,'rel');
