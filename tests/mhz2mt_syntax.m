function [err,data] = test(opt,olddata)

v = mhz2mt;
v = mhz2mt(rand);
v = mhz2mt(rand,rand);
v = mhz2mt(rand,rand(1,6));

err = 0;

data = [];
