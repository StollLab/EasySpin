function [err,data] = test(opt,olddata)

% Multiple spin operators, syntax check

[Sx,Sy,Sz] = sop(1,'x','y','z');
err(1) = 0;

data = [];
