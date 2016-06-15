function [err,data] = test(opt,olddata)

% Syntax test

nu = linspace(9,10,1001);
nu0 = 9.6;
Qu = 2000;
beta = 2;

G = resonatordip(nu,nu0,Qu,beta);
G = resonatordip([],nu0,Qu,beta);

err = false;
data = [];

