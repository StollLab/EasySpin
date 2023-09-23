function ok = test(opt)

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

Vary.g = [0 0.02]; 

FitOpt.Verbosity = 0;
FitOpt.maxGenerations = 30;

FitOpt.Method = [fitAlg ' ' dataMethod];

result = esfit_legacy(spc,@pepper,{Sys,Exp},{Vary},FitOpt);
rmsd = result.rmsd/max(result.fit);

ok = rmsd<0.01;

if opt.Display
  subplot(4,1,[1 2 3]);
  plot(nu,spc,nu,result.fit);
  xlabel('\nu (GHz)');
  legend('sim','fit');
  legend boxoff
  subplot(4,1,4);
  plot(nu,spc-result.fit.');
end
