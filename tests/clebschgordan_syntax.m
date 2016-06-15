function [err,data] = test(opt,olddata)

% Test 1: Input syntax check
%======================================================
a1 = clebschgordan(4,3,2,1,-2,-1);
a2 = clebschgordan([4 3 2],[1 -2 -1]);
a3 = clebschgordan([4 1],[3 -2],[2 -1]);
a4 = -sqrt(7)/6;
err = ~areequal([a1 a2 a3],[a4 a4 a4]);
data = [];
