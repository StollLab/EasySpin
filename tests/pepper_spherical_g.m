function [err,data] = test(opt,olddata)

% Spherical syntax for Sys.g

giso = 2;
gax = 0.1;
grh = 0.03;
lw = 2;

Sys1.g = giso + [-1 -1 2]*gax + [-1 +1 0]*grh;
Sys1.lw = lw;
Sys2.g_ = [giso gax grh];
Sys2.lw = lw;

Exp.mwFreq = 9.8;
Exp.Range = [300 400];

[x,y1] = pepper(Sys1,Exp);
[x,y2] = pepper(Sys2,Exp);

if opt.Display
  plot(x,y1,x,y2);
  legend('cartesian','spherical');
  title(mfilename);
end

err = ~areequal(y1,y2,1e-6,'rel');
data = [];
