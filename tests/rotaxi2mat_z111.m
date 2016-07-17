function [err,data] = test(opt,olddata)

%======================================================
% 120 degree rotations around 111
%======================================================
n = [1;1;1];
rho = 2*pi/3;
R = rotaxi2mat(n,rho);
z = [0;0;1];
y = [0;1;0];
ok = areequal(R*z,y,1e-7);

err = any(~ok);
data = [];
