function [err,data] = test(opt,olddata)

% Isotopologues with isotope-specific rescaling of hyperfine couplings
%-------------------------------------------------------------------------------
A_1H = 10;
A_2H = A_1H*nucgval('2H')/nucgval('1H');

Sys.Nucs = 'H';
Sys.A = A_1H;

Iso = isotopologues(Sys);

ok = Iso(1).A == A_1H && Iso(2).A == A_2H;

err = ~ok;
data = [];
