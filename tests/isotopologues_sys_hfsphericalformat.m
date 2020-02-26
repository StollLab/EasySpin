function ok = test()

% Isotopologues with isotope-specific rescaling of hyperfine couplings
%-------------------------------------------------------------------------------

% isotropic A_
A_1H = 1;
A_2H = A_1H*nucgval('2H')/nucgval('1H');
Sys.Nucs = 'H';
Sys.A_ = A_1H;
Iso = isotopologues(Sys);
ok(1) = all(Iso(1).A_==A_1H) && all(Iso(2).A_==A_2H);

% axial A_
A_1H = [1 2];
A_2H = A_1H*nucgval('2H')/nucgval('1H');
Sys.Nucs = 'H';
Sys.A_ = A_1H;
Iso = isotopologues(Sys);
ok(2) = all(Iso(1).A_==A_1H) && all(Iso(2).A_==A_2H);

% rhombic A_
A_1H = [1 2 3];
A_2H = A_1H*nucgval('2H')/nucgval('1H');
Sys.Nucs = 'H';
Sys.A_ = A_1H;
Iso = isotopologues(Sys);
ok(3) = all(Iso(1).A_==A_1H) && all(Iso(2).A_==A_2H);
