function ok = test()

% Test whether resfields returns intensities that depend on
% photoexcitation.

Triplet.S = 1;
Triplet.D = [900 160];  % MHz
Triplet.lwpp = 1;  % mT
Triplet.Pop = [1 0 1];
Triplet.tdm = 'y';  % transition dipole moment along yMol

Exp.Field = 330;  % mT

Exp.CrystalOrientation = [pi/3 pi/2 pi/8];

Exp.lightMode = '';
[~,Iiso] = resfreqs_matrix(Triplet,Exp);

Exp.lightMode = 'para';
[~,Iperp] = resfreqs_matrix(Triplet,Exp);

Exp.lightMode = 'perp';
[~,Ipara] = resfreqs_matrix(Triplet,Exp);

ok(1) = all(Iiso~=Iperp);
ok(2) = all(Iiso~=Ipara);
ok(3) = all(Iperp~=Ipara);

end