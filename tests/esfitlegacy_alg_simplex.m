function ok = test()

% Assure that running esfit_legacy and pepper using Simplex algorithm is
% successful and yields a good fit.

fitAlg = 'simplex';
dataMethod = 'fcn';

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

spc = pepper(Sys,Exp);
%rng(1)
%data = addnoise(spc,50,'u');

Vary.g = [0.02 0.02];
FitOpt.Verbosity = 0;

FitOpt.Method = [fitAlg ' ' dataMethod];
result = esfit_legacy(spc,@pepper,{Sys,Exp},{Vary},FitOpt);
rmsd = result.rmsd/max(result.fit);

ok = rmsd<1e-10;
