function ok = test()

% Assure saffron doesn't crash when asked to simulate spectra for
% a crystal with space group with multiple sites per unit cell.

Sys.Nucs = '1H';
Sys.A = [5 20];
Sys.lwEndor = 0.5;

Exp2.MolFrame = [pi/3 pi/6 pi/4];
Exp2.SampleFrame = [pi/9 pi/5 0];
Exp2.CrystalSymmetry = 'P212121';

Exp2.Sequence = 'HYSCORE';
Exp2.Field = 330;
Exp2.tau = 0.080;
Exp2.dt = 0.120;
Exp2.nPoints = 256;

Opt.Verbosity = 0;
[x,y] = saffron(Sys,Exp2,Opt);

ok = true;
