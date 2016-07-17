function [err,data] = test(opt,olddata)

a = clight;
b = 299792458;
err = a~=b;

data = [];
