function [err,data] = test(opt,olddata)

%=======================================================
% Test 1: calling syntax
%=======================================================
clear Sys Exp
Sys = struct('S',1/2,'g',[2 3 4],'Nucs','63Cu','A',[56 34 -12]);
Exp = struct('mwFreq',9.5,'Range',[0 600]);
Exp.CrystalOrientation = rand(1,3)*pi;

p = eigfields(Sys,Exp);
[p,i] = eigfields(Sys,Exp);
Opt = [];
p = eigfields(Sys,Exp,Opt);
Opt = struct('Verbosity',0);
p = eigfields(Sys,Exp,Opt);

err = 0;
data = [];
