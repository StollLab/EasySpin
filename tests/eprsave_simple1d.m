function [err,data] = test(opt,olddata)

%======================================================
% Save and load 1D dataset with linear X axis
%======================================================

fn = tempname;

N = 1001;
x = linspace(300,400,N);
y = rand(1,N);
eprsave(fn,x,y);
[xo,yo] = eprload(fn);
delete([fn '.*']);

err = ~areequal(x,xo,1e-5) || ~areequal(y,yo,1e-5);

data = [];

