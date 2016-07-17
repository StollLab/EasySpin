function [err,data] = test(opt,olddata)

% Test 1: Common isotopes
%======================================================
Isotopes = {'1H','2H','13C','14N','15N','33S','59Co','63Cu'};
Spins = [1/2,1,1/2,1,1/2,3/2,7/2,3/2];

for k=1:numel(Isotopes);
  ok = nucspin(Isotopes{k})==Spins(k);
end

err = any(~ok);
data = [];
