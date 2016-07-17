function [err,data] = test(opt,olddata)

% Check whether resfreqs_matrix handles HStrain correctly.

clear Sys Exp
Sys.S = 1/2;
Sys.g = [2 2.05 2.1];
Sys.HStrain = [10 40 100];

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

y = pepper(Sys,Exp);
y = y/max(y);

data.y = y;

if ~isempty(olddata)
  err = ~areequal(y,olddata.y,1e-3);
else
  err = [];
end

