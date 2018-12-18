function [err,data] = test(opt,olddata)

%===============================================================================
% Compare rigid-limit simulation with orientational potential
% using chili and pepper
%===============================================================================

Sys.g = [2 2.2];
Sys.logtcorr = -6;
Sys.lw = [0 10];

Exp.Field = 334;
Exp.mwRange = [9 11];
Exp.CrystalOrientation = [0 0 0];

Opt.LLMK = [100 0 0 0];
Opt.Verbosity = 0;

c200 = 3;

Sys.Potential = [2 0 0 c200];
[B,spc1] = chili(Sys,Exp,Opt);

Exp.Ordering = c200;
[B,spc2] = pepper(Sys,Exp,Opt);

scale = max(spc2);
spc2 = spc2/scale;
spc1 = spc1/scale;

threshold = 0.15;
err = ~areequal(spc1,spc2,threshold);
data = [];

if opt.Display
  subplot(3,1,[1,2])
  plot(B,spc1,B,spc2);
  legend('chili','pepper');
  subplot(3,1,3)
  plot(B,spc1-spc2);
  line(xlim,threshold*[1 1],'Color','r');
  line(xlim,-threshold*[1 1],'Color','r');
end
