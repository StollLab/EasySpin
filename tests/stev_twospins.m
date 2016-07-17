function [err,data] = test(opt,olddata)

% Stevens operators for two-spin system
%======================================================
Spin1 = 5/2;
Spin2 = 3;
Op0 = stev(Spin2,4,3);
Op0 = kron(eye(2*Spin1+1),Op0);
Op1 = stev([Spin1 Spin2],4,3,2);

err = any(abs(Op0(:)-Op1(:))>1e-10);
data = [];
