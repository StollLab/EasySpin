function err = test(opt)

% Test whether rotation via Exp.SampleRotation gives the same results as
% rotation via rotateframe().

Sys.g = [2 2.1 2.3];
Sys.lwpp = 0.2;
Sys.A = [10 50 200];
Sys.AFrame = [0 pi/5 0];
Sys.Nucs = '1H';

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;

nL = [0.2 -0.8 0.5]; % rotation axis, lab-fixed
rho = deg2rad(266); % rotation angle

cori = [0 0 0]; % initial crystal orientation

% Crystal rotation using rotateframe
Exp.CrystalOrientation = rotateframe(cori,nL,rho);
[B,spc2] = pepper(Sys,Exp);

% Crystal rotation using Exp.SampleRotation
Exp.CrystalOrientation = cori;
Exp.SampleRotation = {rho,nL};
[B,spc3] = pepper(Sys,Exp);

err = ~areequal(spc2,spc3,1e-10,'rel');

if opt.Display
  plot(B,spc2,B,spc3)
  axis tight
  legend('rotateframe','Exp.SampleRotation');
end
