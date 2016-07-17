function [err,data] = test(opt,olddata)

% Test 3: a previous bug
%======================================================
Isotopes = {'159Tb'};
Spins = [3/2];

for k = 1:numel(Isotopes)
  ok = nucspin(Isotopes{k})==Spins(k);
end

err = any(~ok);
data = [];
