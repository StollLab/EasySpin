function ok = test()

tau = 0.01;

refExp.Flip = [1 1 2];
refExp.Inc = [0 1 0];
refExp.t = [tau 0 tau];
refExp.dt = 0.012;
refExp.nPoints = 128;
refExp.Phase = [1 1 1];
refExp.tp  = [0 0 0];

Pulse1.Flip = pi/2;
Pulse2.Flip = pi;

newExp.nPoints = refExp.nPoints;

newExp.Sequence = {Pulse1 tau Pulse1 0 Pulse2 tau};
newExp.Dim1 = {'d2' refExp.dt};

Opt.SimulationMode = 'fast';

[~,~,~,processedExp] = runprivate('s_sequencer',newExp,Opt);

ok = all([isequal(refExp,processedExp)]);
