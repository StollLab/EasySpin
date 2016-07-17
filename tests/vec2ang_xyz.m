function [err,data] = test(opt,olddata)

%======================================================
% Principal axes x, y and z
%======================================================
[p,t] = vec2ang([0;0;1]);
ok(1) = areequal([p t],[0 0]);

[p,t] = vec2ang([0;1;0]);
ok(2) = areequal([p t],[pi/2 pi/2]);

[p,t] = vec2ang([1;0;0]);
ok(3) = areequal([p t],[0 pi/2]);

err = any(~ok);
data = [];
