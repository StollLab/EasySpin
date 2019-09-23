function [err,data] = test(opt,olddata)

% Test whether Sys.lw and Sys.lwpp are consistent for
% both Gaussian and Lorentzian broadening

Sys = struct('g',2,'Nucs','1H','lw',[0,0.01],'A',10);
Exp = struct('mwFreq',9.7);

lwG = 0.1; lwL = 0.1;
Sys.lw = [lwG lwL]; Sys.lwpp = [0 0];
[x1,y1] = garlic(Sys,Exp);
Sys.lw = [0 0]; Sys.lwpp = [lwG/sqrt(2*log(2)), lwL/sqrt(3)];
[x2,y2] = garlic(Sys,Exp);

err = ~areequal(y1,y2,1e-10,'rel');

data = [];
