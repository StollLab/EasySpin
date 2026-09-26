function ok = test()

% Resonance fields with higher-order Zeeman terms are the same with sparse and
% full matrices (with transition pre-selection, needed for more than one spin)

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A = 30;
Sys.Ham110 = -sqrt(3)*2*13.9962;  % (lB,lS,l) = (1,1,0), MHz/mT; equivalent to g = 2
Exp.mwFreq = 9.5;
Exp.Range = [300 380];

Opt.Sparse = false;
Bfull = resfields(Sys,Exp,Opt);
Opt.Sparse = true;
Bsparse = resfields(Sys,Exp,Opt);

ok = numel(Bfull)==2 && areequal(Bfull,Bsparse,1e-6,'abs');
