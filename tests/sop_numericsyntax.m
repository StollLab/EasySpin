function ok = test()

% Test numeric syntax, which is used internally

Spins = [3/2 1/2];
threshold = 1e-14;

A = sop(Spins,[1,1]);
B = sop(Spins,'x1');
ok(1) = areequal(A,B,threshold,'abs');

A = sop(Spins,[2,3]);
B = sop(Spins,'z2');
ok(2) = areequal(A,B,threshold,'abs');

A = sop(Spins,[1 4; 2 5]);
B = sop(Spins,'+-');
ok(3) = areequal(A,B,threshold,'abs');
