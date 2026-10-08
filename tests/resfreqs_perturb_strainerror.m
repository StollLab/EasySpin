function ok = test()

% Perturbation-theory solvers do not support strains

Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [10 20 30];
Sys.StrainPars = {'g(1)'};
Sys.StrainFWHM = 0.01;

Exp.mwFreq = 9.5;
Exp.Field = 330;
Exp.Range = [280 360];

fcns = {@resfields_perturb,@resfreqs_perturb,@endorfrq_perturb};
for k = 1:numel(fcns)
  try
    fcns{k}(Sys,Exp);
    ok(k) = false;
  catch ME
    ok(k) = contains(ME.message,'StrainPars');
  end
end
