function ok = test()

tau = 0.01;

RefExp.Field = 350;
RefExp.Flip = [1 1 2];
RefExp.Inc = [0 1 0];
RefExp.t = [tau 0 tau];
RefExp.dt = 0.012;
RefExp.nPoints = 128;


% These are default fields, that are set internally by saffron, but always
% set with s_sequencer and hence need to be added:
RefExp.Phase = [1 1 1];
RefExp.tp  = [0 0 0];

Pulse.Flip = pi/2;
Pulse2.Flip = pi;

NewExp.Sequence = {Pulse tau Pulse 0 Pulse2 tau};
NewExp.Field = RefExp.Field;
NewExp.nPoints = RefExp.nPoints;
NewExp.Dim1 = {'d2' RefExp.dt};

Opt.SimulationMode = 'fast';

[~,~,~,ProcessedExp] = runprivate('s_sequencer',NewExp,Opt);

ProcessedExp = rmfield(ProcessedExp,'Processed');

ok = all([isequal(RefExp,ProcessedExp)]);
