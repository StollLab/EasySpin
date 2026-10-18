function ok = test(opt)

% Check dispersion signals without a convolutional line width, where the phase
% is applied via Hilbert transform: (1) strain-only broadening, compared with
% the phase applied by a narrow Lorentzian convolution, and (2) an isotropic
% spectrum with Gaussian line template, compared with the analytical
% dispersion line shape

Exp.mwFreq = 9.5;
Exp.Range = [315 345];
Exp.nPoints = 4001;
Exp.Harmonic = 0;
Exp.mwPhase = pi/2;

% (1) Strain-only broadening
Sys.g = [2 2.05 2.1];
Sys.lwAnisoG = [50 100 150];
[x,y0] = pepper(Sys,Exp);
Sys.lw = [0 0.03];
[x,y1] = pepper(Sys,Exp);
ok(1) = areequal(y0,y1,1e-2,'rel');

% (2) Isotropic, Gaussian template
Sys = struct('g',2,'lw',[0.5 0]);
x0 = Exp.mwFreq*1e9*planck/bmagn/Sys.g*1e3;  % mT
[x,yDisp] = pepper(Sys,Exp);
Exp.mwPhase = 0;
[x,yAbs] = pepper(Sys,Exp);
yDisp = yDisp/max(yAbs);
yRef = gaussian(x,x0,Sys.lw(1),0,pi/2)/max(gaussian(x,x0,Sys.lw(1)));
ok(2) = areequal(yDisp,yRef,1e-3,'rel');

if opt.Display
  subplot(2,1,1); plot(x,y0,x,y1,'--'); legend('strain only','strain + narrow Lorentzian');
  subplot(2,1,2); plot(x,yDisp,x,yRef,'--'); legend('pepper','analytical');
end

end
