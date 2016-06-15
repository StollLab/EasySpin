function [err,data] = test(opt,olddata)

%====================================================
% syntax check
%====================================================
y = zeros(1,100);
y(1) = -1;
y(2) = 2;
z = 2*y;

y_ = rescale(y,'minmax');
y_ = rescale(y,'maxabs');

y_ = rescale(y,z,'minmax');
y_ = rescale(y,z,'maxabs');
y_ = rescale(y,z,'lsq');
y_ = rescale(y,z,'lsq0');
y_ = rescale(y,z,'lsq1');
y_ = rescale(y,z,'lsq2');
y_ = rescale(y, 'none');

err = 0;
data =[];
