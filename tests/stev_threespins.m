function ok = test()

% Stevens operators for two-spin system

Spin1 = 5/2;
Spin2 = 7/2;
Spin3 = 3/2;
k = 4; q = -3;
Op0 = stev(Spin2,[k,q]);
Op0 = kron(eye(2*Spin1+1),Op0);
Op0 = kron(Op0,eye(2*Spin3+1));
Op1 = stev([Spin1 Spin2 Spin3],[k,q,2]);

ok = areequal(Op0,Op1,1e-10,'abs');
