function ok = test()

% For a predefined experiment (Exp.Sequence is a string/char array),
% s_sequencer should return an empty matrix as 4th output.

refExp.Sequence = '2pESEEM';
refExp.dt = 0.100;
refExp.nPoints = 512;
refExp.Field = 324.9;

Opt.SimulationMode = 'fast';

[~,~,~,processedExp1] = runprivate('s_sequencer',refExp,Opt);

ok(1) = isempty(processedExp1);

% Make sure Options/Additional fields are carried over
refExp2 = refExp;
refExp2.Range = [0 30];
refExp2.mwFreq = 9.4;
refExp2.ExciteWidth = 200;
refExp2.SampleFrame = [0 -pi/2 0];  
refExp2.CrystalSymmetry = 'C2h';

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

[~,~,processedOpt,processedExp2] = runprivate('s_sequencer',refExp2,Opt);

ok(2) = isempty(processedExp2);
ok(3) = isequal(Opt,processedOpt);
