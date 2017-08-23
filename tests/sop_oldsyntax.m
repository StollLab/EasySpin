function [err,data] = test(opt,olddata)

% Old syntax, which is still used internally
%==========================================================

s = [3/2 1/2];

a = sop(s,1,1);
b = sop(s,'xe');

err = ~areequal(a,b,1e-12);

data = [];
