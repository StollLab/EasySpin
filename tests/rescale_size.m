function [err,data] = test(opt,olddata)

%====================================================
% check whether rescale works with vectors of different size
%====================================================
y = zeros(1,100);
y(1) = -1;
y(2) = 2;
z1 = zeros(1,90);
z1(5) = 10;
z2 = zeros(90,1);
z2(5) = 10;

y_ = rescale(y,z1,'minmax');
y_ = rescale(y,z2,'minmax');

err = 0;
data =[];
