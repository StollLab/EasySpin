function [err,data] = test(opt,olddata)

%======================================================
% Save and load 2D dataset with linear X and Y axes
%======================================================

fn = tempname;

n1 = 50;
n2 = 10;
x = 1:n1;
y = 1:n2;
data = rand(n1,n2)+1i*rand(n1,n2);

eprsave(fn,{x,y},data);

[x0,data0] = eprload(fn);

delete([fn '.*']);

thr = 1e-5;
err = ~areequal(x,x0{1},thr,'rel') || ~areequal(y,x0{2},thr,'rel') || ~areequal(data,data0,thr,'rel');

data = [];
