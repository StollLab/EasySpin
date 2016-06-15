function [err,data] = test(opt,olddata)

%======================================================
% Syntax checks
%======================================================
phi = rand;
v = rand(3,1);
rotaxi2mat(1,phi);
rotaxi2mat(2,phi);
rotaxi2mat(3,phi);
rotaxi2mat(v,phi);
rotaxi2mat(v.',phi);
R = rotaxi2mat(3,phi);
R = rotaxi2mat(v,phi);
R = rotaxi2mat(v.',phi);

err = 0;
data = [];
