function [err,data] = test(opt,olddata)

% Test 1: single operator for S=2
%======================================================
Op0 = stev(2,4,2);
a = sqrt(54);
Op1 = [0 0 a 0 0; 0 0 0 -12 0; a 0 0 0 a; 0 -12 0 0 0; 0 0 a 0 0];

err = any(abs(Op0(:)-Op1(:))>1e-10);
data = [];
