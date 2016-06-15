function [err,data] = test(opt,olddata)

%=======================================================
% Test 1: syntax
%=======================================================
Sys = struct('g',2,'Nucs','1H','lw',0.01,'A',10);
Exp = struct('mwFreq',9.7);
Opt = struct('unused',NaN);

y = garlic(Sys,Exp);
[x,y] = garlic(Sys,Exp);
y = garlic(Sys,Exp,Opt);
[x,y] = garlic(Sys,Exp,Opt);

err = 0;
data = [];
