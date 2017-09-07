function [err,data] = test(opt,olddata)

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
rng_(1,'twister')
data = addnoise(spc,50,'u');

Vary.g = [0.02 0.02]; 
Opt = struct;

FitOpt.PrintLevel = 0;

err = false;

rmsd = zeros(1,nMethods);

for iMethod=1:nMethods
  FitOpt.Method = [fitAlg, dataMethod{iMethod}];
  [dummy,dummy,resid] = esfit('pepper',spc,Sys,Vary,Exp,Opt,FitOpt);
  rmsd(iMethod) = sqrt(mean(resid.^2));
end

if any(rmsd>1e-10)
  err = true;
end

data = [];
