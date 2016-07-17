function [err,data] = test(opt,olddata)

Sys.g = rand(3,3)*0.1 + 2;
Sys.lwpp = 0.1;
Exp.mwFreq = 9.5;

[x,y] = garlic(Sys,Exp);

err = 0;
data = [];
