function ok = test()

Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [10 20 40];
Exp.Field = 3500;
Exp.SampleFrame = [0 0 0];

% Rotate sample such that zM, yM and xM end up along zL
% and calculate splitting between ENDOR peaks
rotaxis = [1 1 1];
Exp.SampleRotation = {rotaxis,0};
[nuz,~] = endorfrq_perturb(Sys,Exp);
Exp.SampleRotation = {rotaxis,2*pi/3};
[nuy,~] = endorfrq_perturb(Sys,Exp);
Exp.SampleRotation = {rotaxis,-2*pi/3};
[nux,~] = endorfrq_perturb(Sys,Exp);
hf_splitting = abs([diff(nux) diff(nuy) diff(nuz)]);

% Reference values for ENDOR peak splitting
hf_splitting_ref = Sys.A;

ok = areequal(hf_splitting,hf_splitting_ref,1e-10,'abs');
