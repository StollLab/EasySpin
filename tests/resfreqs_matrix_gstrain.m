function [err,data] = test(opt,olddata)

% Check whether resfreqs_matrix handlex gStrain correctly.


Sys.S = 1/2;
Sys.g = 2;
Sys.gStrain = 0.01;

Exp.Field = 350;
Exp.mwRange = [9.4 10];

[x,y] = pepper(Sys,Exp);
y = y/max(y);

if opt.Display
  plot(x,y,x,olddata.y);
  legend('new','old');
end
data.y = y;

if ~isempty(olddata)
  err = ~areequal(y,olddata.y,1e-3);
else
  err = [];
end

