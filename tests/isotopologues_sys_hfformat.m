function [err,data] = test(opt,olddata)

% Isotopologues with isotope-specific rescaling of hyperfine couplings
%-------------------------------------------------------------------------------

% axial A
A_1H = [1 2];
A_2H = A_1H*nucgval('2H')/nucgval('1H');
Sys.Nucs = 'H';
Sys.A = A_1H;
Iso = isotopologues(Sys);
ok = all(Iso(1).A==A_1H) && all(Iso(2).A==A_2H);

% rhombic A
A_1H = [1 2 3];
A_2H = A_1H*nucgval('2H')/nucgval('1H');
Sys.Nucs = 'H';
Sys.A = A_1H;
Iso = isotopologues(Sys);
ok = all(Iso(1).A==A_1H) && all(Iso(2).A==A_2H);

% full A
A_1H = [1 2 3; 4 5 6; 7 8 9];
A_2H = A_1H*nucgval('2H')/nucgval('1H');
Sys.Nucs = 'H';
Sys.A = A_1H;
Iso = isotopologues(Sys);
ok = all(Iso(1).A(:)==A_1H(:)) && all(Iso(2).A(:)==A_2H(:));

err = ~ok;
data = [];
