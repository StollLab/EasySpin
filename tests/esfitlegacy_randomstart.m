function ok = test(opt)

% Test whether FitOpt.RandomStart works

Sys.g = 2;
Sys.lw = 10;
Exp.mwFreq = 9.5;
Exp.Range = [300 380];
[B,spc] = pepper(Sys,Exp);

Vary.lw = 2;

rng(123);
FitOpt.RandomStart = 1;
FitOpt.Verbosity = 0;

result = esfit_legacy(spc,@pepper,{Sys,Exp},{Vary},FitOpt);

ok = result.rmsd/max(result.fit)<0.01;

if opt.Display
  plot(B,spc,B,result.fit);
  legend('exp','fit');
end
