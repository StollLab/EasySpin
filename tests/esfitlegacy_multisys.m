function ok = test(opt)

% Assure that esfit_legacy runs using multiple systems.

Sys1.g = [2 2.1];
Sys1.lw = 10;
Sys2.g = [1.98 2.03];
Sys2.lw = 5;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

[nu,spc] = pepper({Sys1,Sys2},Exp);
rng(1)
spc = addnoise(spc,80,'u');

Vary1.g = [0.02 0.02]; 
Vary2.g = [0.02 0.02]; 
Opt = struct;
FitOpt.Verbosity = 0;
FitOpt.Method = 'levmar fcn';
result = esfit_legacy(spc,@pepper,{{Sys1,Sys2},Exp,Opt},{{Vary1,Vary2}},FitOpt);

ok = result.rmsd/max(result.fit)<3e-2;

if opt.Display
  plot(nu,spc,nu,result.fit);
  legend('exp','fit');
end

data = [];
