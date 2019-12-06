function [err,data] = test(opt,olddata)

% Assure that running esfit and pepper using Simplex algorithm is
% successful and yields a good fit.

fitAlg = 'simplex';

dataMethod = {'fcn','int','iint','diff','fft'};

nMethods = numel(dataMethod);

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

spc = pepper(Sys,Exp);
%rng_(1,'twister')
%data = addnoise(spc,50,'u');

Vary.g = [0.02 0.02];
FitOpt.PrintLevel = 0;

rmsd = zeros(1,nMethods);

for iMethod = 1:nMethods
  FitOpt.Method = [fitAlg ' ' dataMethod{iMethod}];
  [dummy,dummy,resid] = esfit(@pepper,spc,Sys,Vary,Exp,[],FitOpt);
  rmsd(iMethod) = sqrt(mean(resid.^2));
end

err = any(rmsd>1e-10);

data = [];
