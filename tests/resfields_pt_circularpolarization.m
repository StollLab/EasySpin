function ok = test()

% Test whether linear and circular polarization (both directions) works
% for isotropic g and standard excitation geometry

Exp.Range = [330 350];
Exp.mwFreq = 9.5;
Exp.CrystalOrientation = rand(1,3)*2*pi;

Sys.g = 2;
Exp.mwPolarization = 'linear';
[dum,Int1] = resfields_perturb(Sys,Exp);
Exp.mwPolarization = 'circular+';
[dum,Int2] = resfields_perturb(Sys,Exp);
Exp.mwPolarization = 'circular-';
[dum,Int3] = resfields_perturb(Sys,Exp);

IntPos = [Int1 Int2 Int3];
IntPos = IntPos/max(IntPos);
IntPos_expected = [0.25 1 0];

okPos = all(abs(IntPos-IntPos_expected)<1e-6);

Sys.g = -2;
Exp.mwPolarization = 'linear';
[dum,Int1] = resfields_perturb(Sys,Exp);
Exp.mwPolarization = 'circular+';
[dum,Int2] = resfields_perturb(Sys,Exp);
Exp.mwPolarization = 'circular-';
[dum,Int3] = resfields_perturb(Sys,Exp);

IntNeg = [Int1 Int2 Int3];
IntNeg = IntNeg/max(IntNeg);
IntNeg_expected = [0.25 0 1];

okNeg = all(abs(IntNeg-IntNeg_expected)<1e-6);

ok = okPos && okNeg;
