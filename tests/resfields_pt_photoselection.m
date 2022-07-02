function ok = test()

% Test whether resfields returns intensities that depend on
% photoexcitation.

Triplet.S = 1;
Triplet.D = [900 160];  % MHz
Triplet.lwpp = 1;  % mT
Triplet.Pop = [1 0 1];
Triplet.tdm = 'y';  % transition dipole moment along yMol

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [280 400];  % mT
Exp.Harmonic = 0;

Exp.CrystalOrientation = [pi/3 pi/2 pi/8];

Exp.lightBeam = '';
[~,Iiso] = resfields_perturb(Triplet,Exp);

Exp.lightBeam = 'parallel';
[~,Iperp] = resfields_perturb(Triplet,Exp);

Exp.lightBeam = 'perpendicular';
[~,Ipara] = resfields_perturb(Triplet,Exp);

ok(1) = all(Iiso~=Iperp);
ok(2) = all(Iiso~=Ipara);
ok(3) = all(Iperp~=Ipara);

end
