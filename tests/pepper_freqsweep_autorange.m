function [err,data] = test(opt,olddata)

% Basic frequency sweep with pepper
clear Sys Exp
Sys.g = [2 2.05 2.01];
Sys.lwpp = 10; % MHz

Exp.Field = 340; % mT

[x,y] = pepper(Sys,Exp);

integral = sum(y)*(x(2)-x(1));

if ~isempty(olddata)
  err = max(abs(integral-olddata.integral)/integral)>1e-5;
else
  err = [];
end

if opt.Display
  plot(x,olddata.y,x,y);
end

data.integral = integral;
data.y = y;