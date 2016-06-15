function [err,data] = test(opt,olddata)

% Syntax test
%=======================================================
try
  q = 1000;
  z = rand(1,q);
  y = plegendre(16,-5,z);
  plegendre(16,-5,z);
  plegendre(3,z);
  err = 0;
catch
  err = 1;
end

data = [];
