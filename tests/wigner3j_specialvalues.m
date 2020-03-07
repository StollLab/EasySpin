function ok = test()

% m values don't add up to zero
a = wigner3j([3 4 5; 1 0 0]);
b = 0;
ok(1) = areequal(a,b);

% Two j are 4, one j is zero
a0 = 1/3;
a1 = wigner3j([4 4 0; 0 0 0]);
a2 = wigner3j([4 0 4; 0 0 0]);
a3 = wigner3j([0 4 4; 0 0 0]);
ok(2) = areequal([a1 a2 a3],[a0 a0 a0],1e-10,'abs');

% All j are 4
a = wigner3j([4 4 4; 0 0 0]);
b = 3*sqrt(2/1001);
ok(3) = areequal(a,b,1e-10,'abs');

% All m are zero, sum of j is odd
a = wigner3j([4 5 6; 0 0 0]);
b = 0;
ok(4) = areequal(a,b,1e-10,'abs');

% Triangle condition not satisfied
a = wigner3j([4 5 16; 0 0 0]);
b = 0;
ok(5) = areequal(a,b,1e-10,'abs');
