function ok = test()

% Single operator for S=2

Op0 = stev(2,[4,2]);
a = sqrt(54);
Op1 = [0 0 a 0 0; 0 0 0 -12 0; a 0 0 0 a; 0 -12 0 0 0; 0 0 a 0 0];

ok = areequal(Op0,Op1,1e-10,'abs');
