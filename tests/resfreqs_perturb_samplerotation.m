function ok = test()

Sys.S = 1/2;
Sys.g = [2.0 2.1 2.2];
Exp.Field = 0.35;
Exp.mwRange = [0 500];
Exp.SampleFrame = [0 0 0];

% Rotate sample such that zM, yM and xM end up along zL
% and calculate resonance fields
rotaxis = [1 1 1];
Exp.SampleRotation = {rotaxis,0};
[nuz,~] = resfreqs_perturb(Sys,Exp);
Exp.SampleRotation = {rotaxis,2*pi/3};
[nuy,~] = resfreqs_perturb(Sys,Exp);
Exp.SampleRotation = {rotaxis,-2*pi/3};
[nux,~] = resfreqs_perturb(Sys,Exp);

% Reference values for resonance fields
nuref = unitconvert(Exp.Field,'mT->MHz',Sys.g);

ok = areequal([nux nuy nuz],nuref,1e-10,'abs');
