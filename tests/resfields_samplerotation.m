function ok = test()

Sys.S = 1/2;
Sys.g = [2.0 2.1 2.2];
Exp.mwFreq = 22;
Exp.Range = [0 2000];
Exp.SampleFrame = [0 0 0];

% Rotate sample such that zM, yM and xM end up along zL
% and calculate resonance fields
rotaxis = [1 1 1];
Exp.SampleRotation = {0,rotaxis};
Bz = resfields(Sys,Exp);

Exp.SampleRotation = {2*pi/3,rotaxis};
By = resfields(Sys,Exp);

Exp.SampleRotation = {-2*pi/3,rotaxis};
Bx = resfields(Sys,Exp);

% Reference values for resonance fields
Bx_ref = mhz2mt(Exp.mwFreq*1e3,Sys.g(1));
By_ref = mhz2mt(Exp.mwFreq*1e3,Sys.g(2));
Bz_ref = mhz2mt(Exp.mwFreq*1e3,Sys.g(3));

ok = areequal([Bx By Bz],[Bx_ref By_ref Bz_ref],1e-10,'abs');
