function ok = test()

% Assure that running esfit and pepper using particle swarm optimization
% algorithm is successful and yields a good fit.

rng(20672);

fitAlg = 'swarm';
dataMethod = 'fcn';

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

spc = pepper(Sys,Exp);

Vary.g = [0.02 0.02]; 

FitOpt.PrintLevel = 0;
FitOpt.nParticles = 10;

FitOpt.Method = [fitAlg ' ' dataMethod];
[~,~,residuals] = esfit(@pepper,spc,Sys,Vary,Exp,[],FitOpt);
rmsd = sqrt(mean(residuals.^2));

ok = rmsd<1e-4;
