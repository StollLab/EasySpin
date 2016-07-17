function [err,data] = test(opt,olddata)

%==================================================================
% Test zero-field energies of a and F terms against
% explicit equations
%==================================================================
% Bleaney/Trenam, Proc Roy Soc A, 223, 1-14 (1954)


Sys.S = 5/2;
D = 101;
a = 23;
F = 14.9;
Sys.D = D;
Sys.aF = [a F];
Sys.aFFrame = 3;

e1 = eig(zfield(Sys)); e1 = e1(1:2:6);
e2(1) = D/3 - (a-F)/2 - sqrt((18*D+a-F)^2+80*a^2)/6;
e2(2) = D/3 - (a-F)/2 + sqrt((18*D+a-F)^2+80*a^2)/6;
e2(3) = -2*D/3 + (a-F);
e2 = sort(e2.');

err = any(abs(e2-e1)>1e-8);
data = [];
