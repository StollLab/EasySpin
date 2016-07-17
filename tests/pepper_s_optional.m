function [err,data] = test(opt,olddata)

% The S field in the spin system structure is optional.
% Default is 1/2.

Sys = struct('g',[2 2 2.2],'lw',3);
Exp = struct('mwFreq',9.5,'Range',[280 350]);

y = pepper(Sys,Exp);

err = 0;
data = [];
