function [err,data] = test(opt,olddata)

Par.tcorr = 10*rand()*1e-9;
Par.dt = 0.1e-9;
Par.nSteps = 2000;
Par.nTraj = 1;

[t, R] = stochtraj(Par);

err = false;

data = [];
