function [err,data] = resfields_largehf(opt,olddata)

% very large hfc

A_MHz = 1475.4;
Sys = struct('g',1.999,'Nucs','209Bi','n',[1],'A',A_MHz,'lw',0.5);
Exp = struct('mwFreq',3.8,'nPoints',1e4,'Range',[50 450]);
Exp.CrystalOrientation=[0 0 0];
Opt.Method='matrix';
%pepper(Sys,Exp,Opt);

B=sort(resfields(Sys,Exp));
B0=[158.667; 334.560];

err = any(abs(B-B0)>1e-2);
data = [];
