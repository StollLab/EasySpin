function [err,data] = test(opt,olddata)

% Test whether isotopologue() handles multi-electron system
%-------------------------------------------------------------------------------
Sys.Nucs = '(12,13)C';
Sys.A = 5;
Sys.Abund = [0.9 0.1];

Iso = isotopologues(Sys);

err = ~isempty(Iso(1).Nucs);

data = [];
