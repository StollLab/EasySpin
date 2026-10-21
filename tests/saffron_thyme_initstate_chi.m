function ok = test()

% Powder simulations with Sys.initState: removing electron coherences from
% the initial state must reproduce the explicit average over the third Euler
% angle chi (rotation about the field), which the powder grid omits.
% Tolerance limited by counter-rotating pulse components (nu1/nu0 in the
% simulation frame).

tol = 2e-2;

Pulse.tp = 0.02;
Pulse.Flip = pi/4;
Exp.Sequence = {Pulse 0.05};
Exp.Field = 340;
Exp.mwFreq = 9.53;
Exp.DetWindow = [0 0.05];
Opt.SimulationMode = 'thyme';
Opt.SimFreq = 9.0;

Sys.S = 1/2;
Sys.g = 2.0028;
Sys.Nucs = '1H';
Sys.A = [4 12];
Sys.AFrame = [0.3 0.6 0.2];
Spins = [1/2 1/2];

% Powder simulation with a single orientation (phi,theta), chi = 0
ang = [0.4 0.9];
OptP = Opt;
OptP.GridSymmetry = 'O3';
OptP.GridFrame = [ang 0];

% Explicit average over chi, using crystal orientations
ExpC = Exp;
ExpC.SampleFrame = [0 0 0];
nChi = 3;

% Electron coherence plus nuclear coherence, and nuclear coherence only
% (the latter must be kept, so it tests that only electron coherences are
% removed)
states = {sop(Spins,'xe')+sop(Spins,'zx'), sop(Spins,'zx')};
for s = 1:numel(states)
  Sys.initState = {states{s},'uncoupled'};
  [~,yPowder] = saffron(Sys,Exp,OptP);
  yAvg = 0;
  for k = 1:nChi
    ExpC.MolFrame = -fliplr([ang 2*pi*(k-1)/nChi]);
    [~,y] = saffron(Sys,ExpC,Opt);
    yAvg = yAvg + y/nChi;
  end
  ok(s) = areequal(yPowder,yAvg,tol,'rel');
end
