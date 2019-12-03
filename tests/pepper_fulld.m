function [err,data] = test(opt,olddata)

%=======================================================================
% Different D tensor inputs: ful vs. D only
%=======================================================================
Sys1.S = 1;
Sys2.S = 1;
Sys1.lwpp = 1;
Sys2.lwpp = 1;
d = 400; % MHz
Sys1.D = diag([-1 -1 2]*d/3);
Sys2.D = d;

Exp.mwFreq = 9.5;
Exp.Range = [200 400];

[x,y1] = pepper(Sys1,Exp);
[x,y2] = pepper(Sys2,Exp);

if opt.Display
  plot(x,y1,x,y2,'r');
  legend('Sys.D = d','Sys.D = diag([-1 -1 2]*d/3)');
  legend boxoff
end

err = ~areequal(y1,y2,1e-10,'rel');
data = [];

