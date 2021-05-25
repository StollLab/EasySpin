function ok = test()

j1 = 4;
j2 = 3;
j = 2;
m1 = 1;
m2 = -2;
m = -1;

% Reference value
a = -sqrt(7)/6;

% Use all different input forms
a1 = clebschgordan(j1,j2,2,m1,m2,m);
a2 = clebschgordan([j1 j2 j],[m1 m2 m]);
a3 = clebschgordan([j1 m1],[j2 m2],[j m]);

ok = areequal([a1 a2 a3],[a a a],1e-10,'rel');
