function [err,data] = test(opt,olddata)

% Test 2: some values
%======================================================
a(1) = clebschgordan(11/2,5,9/2,7/2,-3,1/2);
b(1) = -3/2*sqrt(5/143);
err = ~areequal(a,b,1e-10,'rel');

data = [];
