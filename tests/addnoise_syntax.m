function [err,data] = test(opt,olddata)

%======================================================
% Adding noise to data
%======================================================

x = linspace(-1,1,10001);
y = gaussian(x,0,0.3);

yn = addnoise(y,10);
yn = addnoise(y,10,'u');
yn = addnoise(y,10,'n');

data = [];
err = 0;
