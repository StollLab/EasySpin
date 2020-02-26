function ok = test(opt)

% Compare HStrain-broadened spectra for 'matrix' and 'perturb' methods

clear Sys Exp
Sys.g = [2.0 2.1 2.2];
Sys.HStrain = [1 1 1]*200;

Exp.mwFreq = 35;
Exp.Range = [1000 1300];

Opt.Method = 'matrix';
[x,y0] = pepper(Sys,Exp,Opt);
Opt.Method = 'perturb';
[x,y1] = pepper(Sys,Exp,Opt);

if opt.Display
  subplot(2,1,1)
  plot(x,y0,x,y1,'r');
  subplot(2,1,2)
  plot(x,y1-y0);
end

ok = areequal(y0,y1,1e-5,'rel');
