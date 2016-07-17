function [err,data] = test(opt,olddata)

a = eps0;
b = 8.854187817620391e-12; % = 1/mu0/clight^2

err = abs(a-b)/a>1e-10;

data = [];
