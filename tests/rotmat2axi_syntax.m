function [err,data] = test(opt,olddata)

%======================================================
% Syntax tests
%======================================================
R = erot(rand(1,3));

rotmat2axi(R);
[a,b] = rotmat2axi(R);

err = 0;
data = [];
