function [err,data] = test(opt,olddata)

%====================================================
% check whether rescale works with minmax option
%====================================================
N = 100;
y = rand(1,N)-0.9;
yr = rand(1,N)-0.5;
y_ = rescale(y,yr,'maxabs');

thr = 1e-12;
ok = areequal(max(abs(y_)),max(abs(yr)),thr,'abs');
ok = ok && sign(y(1))==sign(y_(1));
err = ~ok;

data = [];