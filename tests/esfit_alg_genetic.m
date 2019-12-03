function [err,data] = test(opt,olddata)

% Assure that running esfit and pepper using genetic algorithm is successful 
% and yields a good fit.

fitAlg = 'genetic';
dataMethod = 'fcn';

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];
Exp.nPoints = 500;

[nu,spc] = pepper(Sys,Exp);
rng_(1,'twister')
%data = addnoise(spc,50,'u');

Vary.g = [0 0.02]; 

FitOpt.PrintLevel = 0;
FitOpt.maxGenerations = 30;

FitOpt.Method = [fitAlg ' ' dataMethod];

[sysFit,spcFit] = esfit(@pepper,spc,Sys,Vary,Exp,[],FitOpt);
rmsd = sqrt(mean((spc-spcFit).^2));
err = any(rmsd/max(spc)>0.01);

if opt.Display
  subplot(4,1,[1 2 3]);
  plot(nu,spc,nu,spcFit);
  xlabel('\nu (GHz)');
  legend('sim','fit');
  legend boxoff
  subplot(4,1,4);
  plot(nu,spc-spcFit);
end

data = [];
