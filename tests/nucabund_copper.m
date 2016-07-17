function [err,data] = test(opt,olddata)

% Test 2: Cu isotopes
%======================================================
w1 = nucabund('63Cu,65Cu');
w2 = [0.6917 0.3083];
err = ~areequal(w1,w2);
data = [];
