function [err,data] = test(opt,olddata)

v = mt2mhz;
v = mt2mhz(rand);
v = mt2mhz(rand,rand);
v = mt2mhz(rand,rand(1,6));

err = 0;

data = [];
