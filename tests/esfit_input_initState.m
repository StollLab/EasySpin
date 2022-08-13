function ok = test(opt)

% Make sure esfit runs with Sys.initState

Sys.S = 1;
Sys.g = 2.0023;
Sys.D = [800 -200];
Sys.lw = 5;

Sys.initState = {[0.2 0.5 0.3],'zerofield'};

Exp.mwFreq = 9.75; % GHz
Exp.Range = [300 400];
Exp.Harmonic = 0;

[B,spc] = pepper(Sys,Exp);
rng(1)
spc = addnoise(spc,100,'n');

Vary.D = [100 100];
Vary.initState{1} = [0 0.2 0.2];

Opt = struct;
FitOpt.Verbosity = 0;
result = esfit(spc,@pepper,{Sys,Exp,Opt},{Vary},FitOpt);

ok = result.rmsd/max(result.fit)<3e-2;

if opt.Display
  plot(B,spc,B,result.fit);
  legend('exp','fit');
end
