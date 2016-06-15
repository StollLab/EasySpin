function [err,data] = test(opt,olddata)

% Test 3: Cl isotopes
%======================================================
w1 = nucabund('35Cl,37Cl');
w2 = [0.7578 0.2422];
err = ~areequal(w1,w2);
data = [];
