function [err,data] = test(opt,olddata)

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 200;
Sys.lw = 10;

Exp.Harmonic = 0;
Exp.Field = 330;
Exp.mwRange = [9.1 9.4];

[x,y] = garlic(Sys,Exp);

data.x = x;
data.y = y;

if ~isempty(olddata)
  ok = areequal(olddata.y,data.y,1e-5,'rel');
  err = ~ok;
else
  err = [];
end



