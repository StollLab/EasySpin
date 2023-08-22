function ok = test()

% Test whether resfields returns intensities that depend on
% photoexcitation.

Triplet.S = 1;
Triplet.D = [900 160];  % MHz
Triplet.lwpp = 1;  % mT
Triplet.initState = {[1 0 1],'zerofield'};
Triplet.tdm = 'y';  % transition dipole moment along yMol

Exp.Field = 330;  % mT

Exp.SampleFrame = [pi/3 pi/2 pi/8];

Exp.lightBeam = '';
[~,Iiso] = resfreqs_matrix(Triplet,Exp);

Exp.lightBeam = 'parallel';
[~,Iperp] = resfreqs_matrix(Triplet,Exp);

Exp.lightBeam = 'perpendicular';
[~,Ipara] = resfreqs_matrix(Triplet,Exp);

ok(1) = all(Iiso~=Iperp);
ok(2) = all(Iiso~=Ipara);
ok(3) = all(Iperp~=Ipara);

end
