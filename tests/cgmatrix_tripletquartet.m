function [ok,data] = test(opt,olddata)

% Regression check of transformation matrix for S1=1 and S2=3/2

U2C = cgmatrix(1,3/2);

data.U2C = U2C;

if ~isempty(olddata)
  ok = areequal(data.U2C,olddata.U2C,1e-10,'abs');
else
  ok = [];
end
