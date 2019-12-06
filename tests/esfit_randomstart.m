function [err,data] = test(opt,olddata)

% Test whether FitOpt.RandomStart works

Sys.g = 2;
Sys.lw = 10;
Exp.mwFreq = 9.5;
Exp.Range = [300 380];
[B,spc] = pepper(Sys,Exp);

Vary.lw = 2;

rng_(123);
FitOpt.RandomStart = 1;
FitOpt.PrintLevel = 0;

[spcfit,dummy,residuals] = esfit(@pepper,spc,Sys,Vary,Exp,[],FitOpt);

rmsd = sqrt(sum(residuals.^2));

err = rmsd>0.01;

if opt.Display
  plot(B,spc,B,spcfit);
  legend('exp','fit');
end

data = [];
