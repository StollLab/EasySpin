function [err,data] = test(opt,olddata)

%=======================================================
% Syntax test
%=======================================================
Sys = struct('g',[2.01 2.02 2.03],'Nucs','14N','A',[10 20 30]);
Sys.tcorr = 1e-8; Sys.lw = 0.1;
Exp = struct('mwFreq',9.7);
Opt = struct('unused',NaN);

y = chili(Sys,Exp);
[x,y] = chili(Sys,Exp);
y = chili(Sys,Exp,Opt);
[x,y] = chili(Sys,Exp,Opt);

err = 0;
data = [];
