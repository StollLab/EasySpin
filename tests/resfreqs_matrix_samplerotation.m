function ok = test()

Sys.S = 1/2;
Sys.g = [2.0 2.1 2.2];
Exp.Field = 0.35;
Exp.mwRange = [0 500];
Exp.SampleFrame = [0 0 0];

% Rotate sample such that zM, yM and xM end up along zL
% and calculate resonance fields
rotaxis = [1 1 1];
Exp.SampleRotation = {0,rotaxis};
[nuz,~] = resfreqs_matrix(Sys,Exp);
Exp.SampleRotation = {2*pi/3,rotaxis};
[nuy,~] = resfreqs_matrix(Sys,Exp);
Exp.SampleRotation = {-2*pi/3,rotaxis};
[nux,~] = resfreqs_matrix(Sys,Exp);

% Reference values for resonance fields
nuref = mt2mhz(Exp.Field,Sys.g);

ok = areequal([nux nuy nuz],nuref,1e-10,'abs');
