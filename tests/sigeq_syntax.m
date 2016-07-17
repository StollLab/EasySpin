function [err,data] = test(opt,olddata)

%======================================================================
% Test calling syntax
%======================================================================
[Sx,Sy,Sz] = sop([1/2 1],'xe','ye','ze');
H = Sx + 5*Sy - 2*Sz;

sigeq(H,10);
sigeq(H,10,'pol');
sigma = sigeq(H,10);
sigma = sigeq(H,10,'pol');

err = 0;
data = [];
