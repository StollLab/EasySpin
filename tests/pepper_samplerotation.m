function ok = test(opt)

% Test whether rotation via Exp.SampleRotation gives the same results as
% rotation via rotateframe().

Sys.g = [2 2.1 2.3];
Sys.lwpp = 0.2;
Sys.A = [10 50 200];
Sys.AFrame = [0 pi/5 0];
Sys.Nucs = '1H';

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;

nL = [0.2 -0.8 0.5];    % rotation axis, lab-fixed
rho = deg2rad(266);     % rotation angle

sampleframe0 = [10 50 -30]*pi/180;   % initial crystal orientation

% Sample rotation using rotateframe()
Exp.SampleFrame = rotateframe(sampleframe0,nL,rho);
[B,spc1] = pepper(Sys,Exp);

% Sample rotation using Exp.SampleRotation
Exp.SampleFrame = sampleframe0;
Exp.SampleRotation = {rho,nL};
[B,spc2] = pepper(Sys,Exp);

ok = areequal(spc1,spc2,1e-10,'rel');

if opt.Display
  plot(B,spc1,B,spc2)
  axis tight
  legend('rotateframe','Exp.SampleRotation');
end
