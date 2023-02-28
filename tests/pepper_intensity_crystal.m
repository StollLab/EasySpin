function ok = test()

% Compare integrated intensity of a powder spectrum of an
% isotropic-g spin system with explicit expressions

Sys.g = [2.1 2.1 2.0];
g = diag(Sys.g);
Sys.lwpp = 1;

Exp.mwFreq = 9.6;
Exp.Range = [320 350];
Exp.Harmonic = 0;
Exp.nPoints = 10000;
phi = 2*pi*rand;
theta = pi*rand;
chi = 2*pi*rand;

Exp.SampleFrame = [-chi -theta -phi];

[x1,y1] = pepper(Sys,Exp);
Int1 = sum(y1)*(x1(2)-x1(1));

R = erot(Exp.SampleFrame);
n0 = R(:,3);
n1 = R(:,1);
g0 = norm(n0.'*g);
u = (n0.'*g)/norm(n0.'*g);
TP = norm(cross(n1.'*g,u))^2;

pre = (bmagn/planck/1e9/2)^2;
dBdE = (planck/bmagn)*1e9/g0;
Int2 = pre*TP*dBdE*(8*pi^2);

ok = areequal(Int1,Int2,1e-3,'abs');
