function ok = test(opt)

% Assure that esfit_legacy runs.

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

[nu,spc] = pepper(Sys,Exp);
rng(1)
spc = addnoise(spc,50,'n');

Vary.g = [0.02 0.02]; 
Opt = struct;
FitOpt.Verbosity = 0;
FitOpt.Method = 'levmar fcn';
result = esfit_legacy(spc,@pepper,{Sys,Exp,Opt},{Vary},FitOpt);

ok = result.rmsd/max(result.fit)<3e-2;

if opt.Display
  plot(nu,spc,nu,result.fit);
  legend('exp','fit');
end
