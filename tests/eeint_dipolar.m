function [err,data] = test(opt,olddata)

%======================================================
% Compare ee and J/dip/dvec for a two-spin system
%======================================================
S = [1/2 1/2];
g = [2 2];

d = 2.3e5; % axial dipolar coupling, MHz
r = 1.87454e4; % rhombic dipolar coupling, MHz

% 1) ee
Sys1.S = S;
Sys1.g = g;
Sys1.ee = d*[1 1 -2] + r*[1 -1 0];
Hee1 = eeint(Sys1);

% 2) dip, 2 values
Sys2.S = S;
Sys2.g = g;
Sys2.dip = [d r];
Hee2 = eeint(Sys2);

% 3) dip, 3 values
Sys3.S = S;
Sys3.g = g;
Sys3.dip = d*[1 1 -2] + r*[1 -1 0];
Hee3 = eeint(Sys3);

% The Hamiltonians should be _numerically_ identical
err = any(Hee1(:)~=Hee2(:)) || any(Hee1(:)~=Hee3(:));

data = [];
