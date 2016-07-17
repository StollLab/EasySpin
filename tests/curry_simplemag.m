function [err,data] = test(opt,olddata)

% Magnetization of a simple coupled spin dimer

Sys.S = [1/2 1/2];
Sys.g = [2 2];
Sys.ee = -2*-4*30e3;

Exp.Temperature = 1;
Exp.Field = linspace(0,17,20)*1e3;

m = curry(Sys,Exp);

data.m = m;

if ~isempty(olddata)
  ok = areequal(data.m,olddata.m,1e-3);
  err = ~ok;
  if opt.Display
    plot(Exp.Field/1e3,olddata.m,Exp.Field/1e3,data.m);
  end
else
  err = [];
end
