function [err,data] = test(opt,olddata)

%======================================================
% Rotations around x, y and z
%======================================================
R = rotaxi2mat(1,pi/2);
ok(1) = areequal(R,[1 0 0; 0 0 1; 0 -1 0],0);

R = rotaxi2mat(2,pi/2);
ok(2) = areequal(R,[0 0 -1; 0 1 0; 1 0 0],0);

R = rotaxi2mat(3,pi/2);
ok(3) = areequal(R,[0 1 0; -1 0 0; 0 0 1],0);

R = rotaxi2mat([1 1 1],2*pi/3);
ok(4) = areequal(R,[0 1 0; 0 0 1; 1 0 0],0);

err = any(~ok);
data = [];
