function ok = test()

% Assure that running esfit and pepper using Levenberg-Marquardt algorithm 
% is successful and yields a good fit.

fitAlg = 'levmar';
dataMethod = 'fcn';

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

spc = pepper(Sys,Exp);
%rng_(1,'twister')
%data = addnoise(spc,50,'u');

Vary.g = [0.02 0.02]; 

FitOpt.PrintLevel = 0;

FitOpt.Method = [fitAlg ' ' dataMethod];
[~,~,resid] = esfit(@pepper,spc,Sys,Vary,Exp,[],FitOpt);
rmsd = sqrt(mean(resid.^2));

ok = rmsd<1e-10;
