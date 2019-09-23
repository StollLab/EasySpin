function [err,data] = test(opt,olddata)

d = degree;
err = ~areequal(d,pi/180,0,'abs');
data = [];
