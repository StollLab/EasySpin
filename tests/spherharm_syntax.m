function [err,data] = test(opt,olddata)

% Syntax test
%=======================================================
try
  q = 1000;
  theta = rand(1,q);
  phi = rand(1,q);
  y = spherharm(5,4,theta,phi);
  err = 0;
catch
  err = 1;
end

data = [];
