function [err,data] = test(opt,olddata)

J = 50*30e3;

Sys.S = [3/2 3/2];

Sys.D = [10 1; 20 6];
Sys.ee = J;
F = spinladder(Sys);

Sys.D = [10 20];
F = spinladder(Sys);

Sys.D = [20 30 40; 10 30 50];
F = spinladder(Sys);

err = 0;

data = [];
