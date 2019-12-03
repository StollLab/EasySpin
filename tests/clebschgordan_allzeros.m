function [err,data] = test(opt,olddata)

a = clebschgordan(0,0,0,0,0,0);
b = 1;
err = ~areequal(a,b,1e-10,'rel');
data = [];
