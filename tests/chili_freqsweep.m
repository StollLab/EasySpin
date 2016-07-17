function [err,data] = test(opt,olddata)

%=======================================================
% Basic frequency sweep
%=======================================================
Sys = struct('g',[2.01 2.02 2.03],'Nucs','14N','A',[10 20 30]);
Sys.tcorr = 1e-8; Sys.lw = 0.1;
Exp = struct('Field',340);
Exp.mwRange = [9.4 9.8];

[x,y] = chili(Sys,Exp);

err = isequal(x(1),Exp.mwRange(1),1e-6) || ...
      isequal(x(end),Exp.mwRange(2),1e-6);
data = [];
