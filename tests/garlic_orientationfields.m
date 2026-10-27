function ok = test()

% Sample orientation fields have no effect on isotropic spectra and are
% ignored without a warning

Sys.g = 2;
Sys.Nucs = '14N';
Sys.A = 40;
Sys.lw = 0.1;
Exp.mwFreq = 9.5;
Exp.Range = [335 343];

[~,y0] = garlic(Sys,Exp);

Exp.SampleFrame = [0 pi/4 0];
Exp.SampleRotation = {'x',pi/4};
Exp.MolFrame = [0 pi/4 0];
Exp.CrystalSymmetry = 'C2h';
lastwarn('');
[~,y1] = garlic(Sys,Exp);

ok(1) = areequal(y1,y0);
ok(2) = isempty(lastwarn);

end
