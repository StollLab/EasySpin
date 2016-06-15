function [err,data] = test(opt,olddata)

% Syntax test
%======================================================
a1 = wigner3j(4,3,2,1,-2,1);
a2 = wigner3j([4 3 2],[1 -2 1]);
a3 = wigner3j([4 1],[3 -2],[2 1]);
a4 = -sqrt(7/5)/6;

err = ~areequal([a1 a2 a3],[a4 a4 a4]);

data = [];
