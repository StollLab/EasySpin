function [err,data] = test(opt,olddata)

% Check whether vec2ang works with inputs with different sizes
%==================================================================
v = rand(3,1);
[p,t] = vec2ang(v);
[p,t] = vec2ang(v.');

v = rand(3,7);
[p,t] = vec2ang(v);
[p,t] = vec2ang(v.');

err = 0;
data = [];
