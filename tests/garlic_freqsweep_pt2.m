function [ok,data] = test(opt,olddata)

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 200;
Sys.lw = 10;

Exp.Harmonic = 0;
Exp.Field = 330;
Exp.mwRange = [9.1 9.4];

Opt.Method = 'perturb2';
[x,y] = garlic(Sys,Exp,Opt);

data.x = x;
data.y = y;

if ~isempty(olddata)
  ok = areequal(olddata.y,data.y,1e-5,'rel');
else
  ok = [];
end



