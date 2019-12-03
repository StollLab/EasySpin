function [err,data] = test(opt,olddata)

x_mT = 1;
g = 2;
v1 = mt2mhz(x_mT,g);
v2 = x_mT*1e-3*g*bmagn/planck/1e6;

err = ~areequal(v1,v2,1e-10,'rel');

data = [];
