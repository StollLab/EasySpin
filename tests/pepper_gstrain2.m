function [err,data] = test(opt,olddata)

% g strain comparison for Method='matrix' and Method='perturb'

Sys = struct('S',1/2,'g',[2 2.1 2.2],'gStrain',[1 2 3]*0.01);
Exp = struct('mwFreq',9.5,'Range',[290 350]);
Opt.Method = 'matrix';
[x,y1] = pepper(Sys,Exp);
Opt.Method = 'perturb';
[x,y2] = pepper(Sys,Exp);

err = ~areequal(y1/max(y1),y2/max(y2),1e-4);
data = [];
