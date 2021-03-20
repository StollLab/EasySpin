function ok = test()

% Assure that running esfit and pepper using grid algorithm is successful 
% and yields a good fit.

fitAlg = 'grid ';

dataMethod = {'fcn','int','iint','diff','fft'};

nMethods = numel(dataMethod);

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

[nu,spc] = pepper(Sys,Exp);
rng(1)

Vary.g = [0.02 0.02]; 
Opt = struct;

FitOpt.PrintLevel = 0;

rmsd = zeros(1,nMethods);

for iMethod = 1:nMethods
  FitOpt.Method = [fitAlg, dataMethod{iMethod}];
  result = esfit(spc,@pepper,{Sys,Exp,Opt},{Vary},FitOpt);
  rmsd(iMethod) = result.rmsd/max(result.fit);
end

ok = rmsd<1e-10;
