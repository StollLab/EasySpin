function [err,data] = test(opt,olddata)

% Test numeric syntax, which is used internally
%==========================================================

Spins = [3/2 1/2];
threshold = 1e-14;

A = sop(Spins,1,1);
B = sop(Spins,'x1');
err(1) = ~areequal(A,B,threshold,'abs');

A = sop(Spins,2,3);
B = sop(Spins,'z2');
err(2) = ~areequal(A,B,threshold,'abs');

A = sop(Spins,[1 2],[4 5]);
B = sop(Spins,'+-');
err(3) = ~areequal(A,B,threshold,'abs');

data = [];
