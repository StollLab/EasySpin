function ok = test()

% Compare integrated intensity of a powder spectrum of an
% isotropic-g spin system with explicit expressions

Sys.g = [2];
Sys.lwpp = 0.5;

Exp.mwFreq = 9.6;
Exp.Range = [341 345];
Exp.Harmonic = 0;

Opt.Method = 'perturb2';
[x,y] = pepper(Sys,Exp,Opt);
Int = sum(y)*(x(2)-x(1));

Rate = (bmagn/planck/1e9/2)^2 * Sys.g^2;
dBdE = (planck/bmagn)*1e9/mean(Sys.g);
Int2 = Rate*dBdE;
Int2 = Int2*(8*pi^2);

ok = areequal(Int,Int2,1e-4,'abs');
