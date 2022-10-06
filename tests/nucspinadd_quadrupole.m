function ok = test()

% Check identity between 2-element Q and 3-element Q

A = 1;
eeqQ = 1;
eta = 0.1;
I = 1;
Q = eeqQ/(4*I*(2*I-1)) * [-1+eta, -1-eta, 2];

Sys = struct('S',1/2,'g',[2 2 2.2]);
Sys1 = nucspinadd(Sys,'14N',A,[],[eeqQ eta],[]);
Sys2 = nucspinadd(Sys,'14N',A,[],Q,[]);

H1 = ham_nq(Sys1);
H2 = ham_nq(Sys2);

ok = all(H1(:)==H2(:));
