function ok = test()

Sys.S = 1/2;
Sys.g = [2.0 2.1 2.2];
Exp.mwFreq = 22;
Exp.Range = [0 2000];
Exp.CrystalOrientation = [0 0 0];

% Rotate sample such that zM, yM and xM end up along zL
% and calculate resonance fields
rotaxis = [1 1 1];
Exp.SampleRotation = {0,rotaxis};
[Bz,~] = resfields(Sys,Exp);
Exp.SampleRotation = {2*pi/3,rotaxis};
[By,~] = resfields(Sys,Exp);
Exp.SampleRotation = {-2*pi/3,rotaxis};
[Bx,~] = resfields(Sys,Exp);

% Reference values for resonance fields
Bref = mhz2mt(Exp.mwFreq*1e3,Sys.g);

ok = areequal([Bx By Bz],Bref,1e-10,'abs');
