function [err,data] = test(opt,olddata)

v = rand(3,10);
[p,t] = vec2ang(v);
a = vec2ang(v);

err = any(a~=[p;t]);
data = [];
