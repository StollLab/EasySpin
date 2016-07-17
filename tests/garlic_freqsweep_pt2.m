function [err,data] = test(opt,olddata)

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 200;
Sys.lw = 10;

Exp.Harmonic = 0;
Exp.Field = 330;
Exp.mwRange = [9.1 9.4];

Opt.Method = 'perturb2';
[x,y] = garlic(Sys,Exp,Opt);
ymax = max(y);

data.x = x;
data.y = y;

if ~isempty(olddata)
  ok = areequal(olddata.y,data.y,ymax*1e-5);
  err = ~ok;
else
  err = [];
end



