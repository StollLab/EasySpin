function [err,data] = test(opt,olddata)

% Test intensity of garlic w tcorr against garlic with lw
%-------------------------------------------------------------
g = [2.0 2.1 2.2];
tcorr = 0.01e-9;

Sys2.g = g;
Sys2.tcorr = tcorr;

lwL = 0.362275;
Sys1.g = g;
Sys1.lw = [0 lwL];

Exp.mwFreq = 9.5;
Exp.Range = [320 328];
Exp.Harmonic = 0;

[x,y1] = garlic(Sys1,Exp);
[x,y2] = garlic(Sys2,Exp);

dx = x(2)-x(1);
int1 = sum(y1)*dx;
int2 = sum(y2)*dx;

err = ~areequal(int1,int2,1e-3*int1);
data = [];

if opt.Display
  [int1 int2]
  plot(x,y1,x,y2);
  legend('garlic lw','garlic tcorr');
end

