function [err,data] = test(opt,olddata)

% Interface test
%-------------------------------------------------------
y = rand(1,500);
m = 5;
yy = datasmooth(y,m);
yy = datasmooth(y,m,'savgol');
yy = datasmooth(y,m,'binom');
yy = datasmooth(y,m,'flat');
yy = datasmooth(y,m,'savgol',3);
yy = datasmooth(y,m,'savgol',3,1);

err = 0;
data = [];
