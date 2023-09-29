function ok = test()

% Assure that running esfit_legacy and pepper using grid algorithm is successful 
% and yields the correct fit.

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

spc = pepper(Sys,Exp);

Vary.g = [0.02 0.02]; 
Opt = struct;

FitOpt.Verbosity = 0;
FitOpt.GridSize = 3;
FitOpt.RandomizeGrid = false;

dataMethod = {'fcn','int','dint','diff','fft'};
nMethods = numel(dataMethod);
rmsd = zeros(1,nMethods);
for iMethod = 1:nMethods
  FitOpt.Method = ['grid ', dataMethod{iMethod}];
  result = esfit_legacy(spc,@pepper,{Sys,Exp,Opt},{Vary},FitOpt);
  rmsd(iMethod) = result.rmsd/max(result.fit);
end

ok = rmsd<1e-10;
