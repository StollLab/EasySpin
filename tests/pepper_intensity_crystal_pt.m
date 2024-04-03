function ok = test()

% Compare line intensity (perturbation theory) of a single-orientation spectrum of an
% spin system with anisotrpic g with an explicit analytical expression.

Sys.g = [2.1 2.1 2.0];
g = diag(Sys.g);
Sys.lwpp = 1;

Exp.mwFreq = 9.6;
Exp.Range = [320 350];
Exp.Harmonic = 0;
Exp.nPoints = 10000;

alpha = 2*pi*rand;
beta = pi*rand;
gamma = 2*pi*rand;
Exp.MolFrame = [alpha beta gamma];

Opt.Method = 'perturb2';
[x1,y1] = pepper(Sys,Exp,Opt);
Int1 = sum(y1)*(x1(2)-x1(1));

R = erot(Exp.MolFrame);
n0 = R(:,3);  % B0 field direction (lab z) in molecular frame
n1 = R(:,1);  % B1 field direction (lab x) in molecular frame
g0 = norm(n0.'*g);
u = (n0.'*g)/norm(n0.'*g);
TP = norm(cross(n1.'*g,u))^2;

pre = (bmagn/planck/1e9/2)^2;
dBdE = (planck/bmagn)*1e9/g0;
Int2 = pre*TP*dBdE*(8*pi^2);

ok = areequal(Int1,Int2,1e-3,'abs');
