function [err,data] = test(opt,olddata)

xMHz = 1;
g = 2;
v1 = mhz2mt(xMHz,g);
v2 = xMHz*1e6*planck/g/bmagn/1e-3;

err = ~areequal(v1,v2,1e-10,'rel');

data = [];
