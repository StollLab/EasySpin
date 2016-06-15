function [err,data] = test(opt,olddata)

a = avogadro;
b = 6.022140857e23;
if abs(a-b)/b > 1e-10;
  err = 1;
else
  err = 0;
end

data = [];
