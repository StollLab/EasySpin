function [err,data] = test(opt,olddata)

% Test 1: N isotopes
%======================================================
w1 = nucabund('14N,15N');
w2 = [0.99632 0.00368];
err = ~areequal(w1,w2,1e-10,'rel');
data = [];
