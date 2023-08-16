function ok = test(opt)

% Make sure pepper returns orientations separately with Opt.Output = 'separate'
% in a crystal simulation

Sys.g = [2 2.1 2.2];
Sys.lwpp = 0.5;

ori1 = [0 0 0];
ori2 = [10 45 30]*pi/180;

Exp.mwFreq = 9.5;
Exp.Range = [300 400];
Exp.MolFrame = [0 0 0];

Opt.Output = 'separate';
Exp.SampleFrame = [ori1; ori2];
[B,spc3] = pepper(Sys,Exp,Opt);

Opt.Output = 'summed';
Exp.SampleFrame = ori1;
spc1 = pepper(Sys,Exp,Opt);
Exp.SampleFrame = ori2;
spc2 = pepper(Sys,Exp,Opt);

if opt.Display
  subplot(2,1,1)
  plot(B,spc3)
  title('Opt.Output=''separate''');
  xlabel('magnetic field (mT)');
  subplot(2,1,2)
  plot(B,[spc1;spc2]/2)
  title('two separate sims');
  xlabel('magnetic field (mT)');
end

ok(1) = size(spc3,1)==2;
ok(2) = areequal([spc1;spc2]/2,spc3,1e-10,'abs');
