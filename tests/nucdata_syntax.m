function [err,data] = test(opt,olddata)

% Test 1: Syntax
%======================================================
Isotopes = '112Sn';
[s] = nucdata(Isotopes);
[s,gn] = nucdata(Isotopes);
[s,gn,qm] = nucdata(Isotopes);
[s,gn,qm,ab] = nucdata(Isotopes);

err = 0;
data = [];
