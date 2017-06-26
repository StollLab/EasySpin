function [err,data] = test(opt,olddata)

% Assure that running esfit and pepper using Monte Carlo algorithm is 
% successful and yields a good fit using bounds determined by seeded rng 
% results.

fitAlg = 'montecarlo ';

dataMethod = {'fcn','int','iint','diff','fft'};

nMethods = numel(dataMethod);

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

[nu,spc] = pepper(Sys,Exp);
rng_(1,'twister')
data = addnoise(spc,50,'u');

Vary.g = [0.01 0.01]; 
Opt = struct;

FitOpt.PrintLevel = 0;
FitOpt.nTrials = 200;

err = false;

rmsd = zeros(1,nMethods);

for iMethod=1:nMethods
  FitOpt.Method = [fitAlg, dataMethod{iMethod}];
  [dummy,dummy,resid] = esfit('pepper',spc,Sys,Vary,Exp,Opt,FitOpt);
  rmsd(iMethod) = sqrt(mean(resid.^2));
end

errors = [0.01, 1, 100, 1, 0.001];  % determined empirically based on 
                                    % results after seeding

if any(rmsd>errors)
  err = true;
end

data = [];
