function [err,data] = test(opt,olddata)

%======================================================
% Check whether 111-direction gives correct angles
%======================================================
[p,t] = vec2ang([1 1 1]);

err = ~areequal([p t],[pi/4 acos(sqrt(1/3))],1e-12,'rel');
data = [];
