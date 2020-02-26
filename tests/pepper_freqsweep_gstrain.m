function [ok,data] = test(opt,olddata)

% Frequency sweep with pepper, g strain
clear Sys Exp
Sys.g = [2.05 2 1.95];
Sys.gStrain = [0.02 0.01 0.003];

Exp.Field = 340; % mT
Exp.mwRange = [9 10]; % GHz

[x,y] = pepper(Sys,Exp);

data.y = y;

if ~isempty(olddata)
  ok = areequal(y/max(y),olddata.y/max(olddata.y),1e-4,'abs');
else
  ok = [];
end
