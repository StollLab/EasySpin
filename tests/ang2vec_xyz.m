function [err,data] = test(opt,olddata)

% Principal axes

x = ang2vec(0,pi/2);
x0 = [1;0;0];
y = ang2vec(pi/2,pi/2);
y0 = [0;1;0];
z = ang2vec(0,0);
z0 = [0;0;1];

err = ~areequal([x y z],[x0 y0 z0]);

data = [];
