function [err,data] = test(opt,olddata)

a = molgas;
b = avogadro*boltzm;
err = ~areequal(a,b);
data = [];
