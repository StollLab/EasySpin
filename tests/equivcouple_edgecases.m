function ok = test()

% n=0, n=1 and integer spins I>=2

% n=0: a single spin-0
[F,N] = equivcouple(1/2,0);
ok(1) = isequal(F,0) && isequal(N,1);

% n=1: the spin itself, with no zero-abundance entries
[F,N] = equivcouple(1,1);
ok(2) = isequal(F,1) && isequal(N,1);
[F,N] = equivcouple(5/2,1);
ok(3) = isequal(F,5/2) && isequal(N,1);

% I=0
[F,N] = equivcouple(0,3);
ok(4) = isequal(F,0) && isequal(N,1);

% I=2, n=3
I = 2;
n = 3;
[F,N] = equivcouple(I,n);
ok(5) = isequal(F,6:-1:0) && isequal(N,[1 2 3 4 5 3 1]);
ok(6) = sum((2*F+1).*N)==(2*I+1)^n;

end
