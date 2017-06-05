function [err,data] = test(opt,olddata)

Sys.Nucs = '51V,1H';
Sys.A_ = [100 30; 20 2];

H1 = sham(Sys,[0 0 350]);
Sys = nucspinrmv(Sys,2);
H2 = sham(Sys,[0 0 350]);

err = length(H2)~=16;
err = err || size(Sys.A_,1)~=1 || size(Sys.A_,2)~=2;

data = [];
