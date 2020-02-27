function ok = test(opt)

% Assure that esfit runs.

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

[nu,spc] = pepper(Sys,Exp);
rng(1)
spc = addnoise(spc,50,'u');

Vary.g = [0.02 0.02]; 
Opt = struct;
FitOpt.PrintLevel = 0;
FitOpt.Method = 'levmar fcn';
[~,spcfit,resid] = esfit(@pepper,spc,Sys,Vary,Exp,Opt,FitOpt);

rmsd = sqrt(mean(resid.^2));

ok = rmsd<3e-2;

if opt.Display
  plot(nu,spc,nu,spcfit);
  legend('exp','fit');
end
