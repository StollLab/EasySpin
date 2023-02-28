function ok = test()

% This definition of Exp should return the experiment structure for the old saffron
% syntax, named RefExp

% The default call:
RefExp.Sequence = '2pESEEM';
RefExp.dt = 0.100;
RefExp.nPoints = 512;
RefExp.Field = 324.9;

Opt.SimulationMode = 'fast';

[~,~,~,ProcessedExp1] = runprivate('s_sequencer',RefExp,Opt);

ProcessedExp1 = rmfield(ProcessedExp1,'Processed');

% Make sure Options/Additional fields are carried over
RefExp2 = RefExp;
RefExp2.Range = [0 30];
RefExp2.mwFreq = 9.4;
RefExp2.ExciteWidth = 200;
RefExp2.SampleFrame = [0 -pi/2 0];  
RefExp2.CrystalSymmetry = 'C2h';

Opt.GridSize = 20;
Opt.TimeDomain = 0;
Opt.Expand = 2;
Opt.ProductRule = 1;
Opt.nOffsets = 10;
Opt.EndorMethod = 2;
Opt.Window = 'bar';
Opt.lwOffset = 100;
Opt.logplot = 1;
Opt.ZeroFillFactor = 4;

[~, ~,ProcessedOpt,ProcessedExp2] = runprivate('s_sequencer',RefExp2,Opt);

ProcessedExp2 = rmfield(ProcessedExp2,'Processed');

ok = all([isequal(RefExp,ProcessedExp1) isequal(RefExp2,ProcessedExp2) isequal(Opt,ProcessedOpt)]);
