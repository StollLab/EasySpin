function [err,data] = test(opt,olddata)

% Test 2: result
%======================================================
a = rand(1,3);
R1 = erot(a);
b = eulang(R1);
err = ~areequal(a,b);

data = [];
