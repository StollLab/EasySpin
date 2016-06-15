function [err,data] = test(opt,olddata)

% Input and output syntax check
%======================================================
S1 = 2;
S2 = 3/2;
N = hsdim([S1 S2]);
C = cgmatrix(S1,S2);
err = numel(C)~=N^2;
data = [];
