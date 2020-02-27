function ok = test()

% Comparison of two output modes: full matrix vs. rows

rng(554);
angles = rand(1,3);
R1 = erot(angles);
[x,y,z] = erot(angles,'cols');
R2 = [x y z];

ok = areequal(R1,R2);
