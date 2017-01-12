function [err,data] = test(opt,olddata)

% Assure that esfit runs.

Sys.g = [2 2.1];
Sys.lw = 10;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

[nu,spc] = pepper(Sys,Exp);
spc = addnoise(spc,50,'u');

Vary.g = [0.02 0.02]; 
Opt = struct;
FitOpt.PrintLevel = 0;
FitOpt.Method = 'levmar fcn';
[bestsys,spcfit] = esfit('pepper',spc,Sys,Vary,Exp,Opt,FitOpt);

if opt.Display
  plot(nu,spc,nu,spcfit);
  legend('exp','fit');
end

data = [];
err = false;
