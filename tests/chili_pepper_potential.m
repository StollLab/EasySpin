function ok = test(opt)

% Compare rigid-limit simulation with orientational potential
% using chili and pepper

Sys.g = [2 2.2];
Sys.logtcorr = -6;
Sys.lw = [0 10];

Exp.Field = 334;
Exp.mwRange = [9 11];

Opt.LLMK = [100 0 0 0];
Opt.Verbosity = 0;
Opt.GridSymmetry = 'Dinfh';

c200 = 3;

Exp.Ordering = c200;
[B,spc2] = pepper(Sys,Exp,Opt);

Exp.SampleFrame = [0 0 0];
Sys.Potential = [2 0 0 c200];
[B,spc1] = chili(Sys,Exp,Opt);


spc2 = spc2/max(spc2);
spc1 = spc1/max(spc1);

threshold = 0.15;
ok = areequal(spc1,spc2,threshold,'abs');

if opt.Display
  subplot(3,1,[1,2])
  plot(B,spc1,B,spc2);
  legend('chili','pepper');
  subplot(3,1,3)
  plot(B,spc1-spc2);
  yline(threshold,'Color','r','LineStyle','--');
  yline(-threshold,'Color','r','LineStyle','--');
end
