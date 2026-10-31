function ok = test()

% Predefined sequences in thyme mode

Sys.S = 1/2;
Sys.g = 2;

Exp.Field = 350;
Exp.mwFreq = 9.8;
Exp.tau = 0.1;
Exp.dt = 0.02;

Opt.GridSize = 7;
Opt.SimulationMode = 'thyme';

% HYSCORE with a single number of points: used for both dimensions
Exp.Sequence = 'HYSCORE';
Exp.t1 = 0.1;
Exp.t2 = 0.1;
Exp.nPoints = 4;
[x,y] = saffron(Sys,Exp,Opt);
ok(1) = isequal(size(y),[4 4]) && numel(x{1})==4 && numel(x{2})==4;

% Two-pulse ESEEM with Exp.DetSequence (replaced by Exp.DetWindow)
Exp = rmfield(Exp,{'t1','t2'});
Exp.Sequence = '2pESEEM';
Exp.DetSequence = [0 1];
[~,y] = saffron(Sys,Exp,Opt);
ok(2) = numel(y)==4;
