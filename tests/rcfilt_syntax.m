function [err,data] = test(opt,olddata)

y = rcfilt(rand(1,101),1,12);
err = 0;

data = [];
