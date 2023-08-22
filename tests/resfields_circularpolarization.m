function ok = test()

% Test whether linear and circular polarization (both directions) works
% for isotropic g and standard excitation geometry

Exp.Range = [330 350];
Exp.mwFreq = 9.5;
Exp.SampleFrame = rand(1,3)*2*pi;

Sys.g = 2;
Exp.mwMode = {pi/2 pi/2};
[dum,Int1] = resfields(Sys,Exp);
Exp.mwMode = {0 'circular+'};
[dum,Int2] = resfields(Sys,Exp);
Exp.mwMode = {0 'circular-'};
[dum,Int3] = resfields(Sys,Exp);

IntPos = [Int1 Int2 Int3];
IntPos = IntPos/max(IntPos);
IntPos_expected = [0.25 1 0];

ok(1) = all(abs(IntPos-IntPos_expected)<1e-6);

Sys.g = -2;
Exp.mwMode = {pi/2 pi/2};
[dum,Int1] = resfields(Sys,Exp);
Exp.mwMode = {0 'circular+'};
[dum,Int2] = resfields(Sys,Exp);
Exp.mwMode = {0 'circular-'};
[dum,Int3] = resfields(Sys,Exp);

IntNeg = [Int1 Int2 Int3];
IntNeg = IntNeg/max(IntNeg);
IntNeg_expected = [0.25 0 1];

ok(2) = all(abs(IntNeg-IntNeg_expected)<1e-6);
