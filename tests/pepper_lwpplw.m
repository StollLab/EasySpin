function [err,data] = test(opt,olddata)

% Test whether Sys.lw and Sys.lwpp are consistent for
% both Gaussian and Lorentzian broadening.

Sys = struct('g',2,'Nucs','1H','lw',[0,0.01],'A',10);
Exp = struct('mwFreq',9.7,'CenterSweep',[346.5 1]);

lwG = 0.1; lwL = 0.1;
Sys.lw = [lwG lwL]; Sys.lwpp = [0 0];
[x1,y1] = pepper(Sys,Exp);
Sys.lw = [0 0]; Sys.lwpp = [lwG/sqrt(2*log(2)), lwL/sqrt(3)];
[x2,y2] = pepper(Sys,Exp);

%plot(x1,y1,'r',x2,y2,'g');
%legend('lw','lwpp');

err = ~areequal(y1,y2);

data = [];
