function [err,data] = test(opt,olddata)

%============================================================================
% Check whether single-matrix-element call to wignerd() works with an array
% of angles.
%============================================================================

rng_(62323,'twister');
n1 = 3;
n2 = 10;

alpha = rand(n1,n2)*2*pi;
beta = rand(n1,n2)*pi;
gamma = rand(n1,n2)*2*pi;

J = 3; m1 = -2; m2 = 1;

D = wignerd([J m1 m2],alpha,beta,gamma);

err = size(D,1)~=n1 || size(D,2)~=n2;

data = [];