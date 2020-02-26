function ok = test()

% Syntax tests

R = erot(rand(1,3));

rotmat2axi(R);
[a,b] = rotmat2axi(R);

ok = true;
