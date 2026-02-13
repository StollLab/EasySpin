function ok = test()

% Test whether g strain linewidths are correct along principal axes

Sys.g = [2 2 2];
Sys.gStrain = [0.01 0.02 0.04];
Exp.mwFreq = 9.5;
Exp.Range = [300 380];

Exp.SampleFrame = [0 0 0];
[Bz,~,Wz] = resfields_perturb(Sys,Exp);
Exp.SampleFrame = [0 pi/2 0];
[Bx,~,Wx] = resfields_perturb(Sys,Exp);
Exp.SampleFrame = [0 pi/2 pi/2];
[By,~,Wy] = resfields_perturb(Sys,Exp);

ok = areequal([Wx/Bx Wy/By Wz/Bz],Sys.gStrain./Sys.g,1e-10,'rel');

end