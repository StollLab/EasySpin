function [err,data] = test(opt,olddata)

%======================================================
% wignerd syntax
%======================================================

J = 15/2;
angles = rand(1,3)*pi;
wignerd(J,angles);
wignerd(J,angles,'+');
wignerd(J,angles,'-');

beta = rand*pi;
wignerd(J,beta);
wignerd(J,beta,'+');
wignerd(J,beta,'-');

err = 0;
data = [];