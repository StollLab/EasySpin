function [err,data] = test(opt,olddata)

Sys.Nucs = '1H,15N,13C';
A1 = rand(3)*10;
A2 = rand(3)*3;
A3 = rand(3)*6;
Sys.A = [A1;A2;A3];

nucspinrmv(Sys,[]);
nucspinrmv(Sys,1);
nucspinrmv(Sys,2);
nucspinrmv(Sys,1);
nucspinrmv(Sys,[1 2]);
nucspinrmv(Sys,[1 3]);
nucspinrmv(Sys,[2 3]);
nucspinrmv(Sys,[1 2 3]);

err = 0;
data = [];
