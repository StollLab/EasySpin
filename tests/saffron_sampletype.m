function ok = test() %#ok<FNDEF>

% Powder or crystal, as in pepper: a crystal is simulated if
% Exp.CrystalSymmetry or Exp.MolFrame is given, otherwise a powder, for
% which Exp.SampleFrame is ignored. (Averaging over crystal sites is
% tested in saffron_oriweights.)

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz
Sys.AFrame = [0 pi/3 0];

Exp.Field = 350;  % mT
Exp.Sequence = '2pESEEM';
Exp.tau = 0.1;  % µs
Exp.dt = 0.02;  % µs
Exp.nPoints = 64;

Opt.GridSize = 19;

% Powder: Exp.SampleFrame is ignored
[~,yPowder] = saffron(Sys,Exp,Opt);
E = Exp;
E.SampleFrame = [0 pi/2 pi/4];
[~,ySampleFrame,info] = saffron(Sys,E,Opt);
ok(1) = areequal(ySampleFrame,yPowder,1e-12,'abs');
ok(2) = strcmp(info.Method,'subspace');

% Single orientation: Exp.MolFrame or Exp.CrystalSymmetry
E = Exp;
E.MolFrame = [0 0 0];
[~,yMolFrame] = saffron(Sys,E,Opt);
E = Exp;
E.CrystalSymmetry = 'P1';
[~,ySymmetry] = saffron(Sys,E,Opt);
ok(3) = areequal(yMolFrame,ySymmetry,1e-12,'abs');
ok(4) = ~areequal(yMolFrame,yPowder,1e-3,'abs');

% The crystal orientation (Exp.SampleFrame) matters for a crystal
E = Exp;
E.MolFrame = [0 0 0];
E.SampleFrame = [0 pi/2 0];
[~,yTilted] = saffron(Sys,E,Opt);
ok(5) = ~areequal(yTilted,yMolFrame,1e-3,'abs');

end
