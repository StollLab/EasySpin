function [err,data] = test(opt,olddata)

%==================================================================
% Test Bkq terms
%==================================================================
S = 7/2;
B43 = rand;

Sys.S = S;
Sys.B4 = [0 B43 0 0, 0, 0 0 0 0];

Op1 = zfield(Sys);
Op2 = B43*stev(S,4,3);

err = any(abs(Op1(:)-Op2(:))>1e-10);
data = [];
