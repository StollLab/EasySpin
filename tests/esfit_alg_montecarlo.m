function ok = test()

% Assure that running esfit and pepper using Monte Carlo algorithm is 
% successful and yields a good fit using bounds determined by seeded rng 
% results.

fitAlg = 'montecarlo';
dataMethod = 'fcn';

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

spc = pepper(Sys,Exp);
rng(100)
%data = addnoise(spc,50,'u');

Vary.g = [0.01 0.01]; 
Opt = struct;

FitOpt.PrintLevel = 0;
FitOpt.nTrials = 300;

FitOpt.Method = [fitAlg ' ' dataMethod];
result = esfit(spc,@pepper,{Sys,Exp,Opt},{Vary},FitOpt);
rmsd = result.rmsd/max(result.fit);

ok = rmsd<0.005;
