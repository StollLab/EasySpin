function ok = test()

% Fit a strain width with esfit: Sys.StrainPars (cell array) is passed
% through, and Sys.StrainFWHM is fitted.

Sys.g = [2 2.1 2.2];
Sys.StrainPars = {'g(1)','g(3)'};
Sys.StrainFWHM = [0.02 0.03];
Exp.mwFreq = 9.5;
Exp.Range = [295 350];
Exp.nPoints = 500;
Opt.GridSize = 10;

[~,spc] = pepper(Sys,Exp,Opt);

Sys0 = Sys;
Sys0.StrainFWHM = [0.015 0.035];
Vary.StrainFWHM = [0.01 0.01];

FitOpt.Method = 'levmar fcn';
FitOpt.Verbosity = 0;
result = esfit(spc,@pepper,{Sys0,Exp,Opt},{Vary},FitOpt);

ok(1) = areequal(result.pfit,Sys.StrainFWHM(:),1e-3,'abs');
ok(2) = isequal(result.argsfit{1}.StrainPars,Sys.StrainPars);
