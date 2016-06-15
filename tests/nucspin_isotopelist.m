function [err,data] = test(opt,olddata)

% Test 2: Isotope list
%======================================================
Isotopes = '1H,1H,2H,15N';
Spins = [1/2 1/2 1 1/2];
s = nucspin(Isotopes);
ok = s==Spins;
err = any(~ok);
data = [];
