function [ok,data] = test(opt,olddata)

% Basic frequency sweep with pepper
clear Sys
Sys.g = [2 2.05 2.01];
Sys.lwpp = 10; % MHz

Exp.Field = 340; % mT
Exp.mwRange = [9 10]; % GHz

[x,y] = pepper(Sys,Exp);

data.y = y;

if ~isempty(olddata)
  ok = areequal(y/max(y),olddata.y/max(olddata.y),1e-4,'abs');
else
  ok = [];
end
