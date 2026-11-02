function ok = test()

% Predefined experiments have a fixed detection phase: a nonzero
% Exp.DetPhase is rejected in both simulation modes

Sys.Nucs = '1H';
Sys.A = 5;
Exp.Field = 340;
Exp.Sequence = '2pESEEM';
Exp.dt = 0.01;
Exp.DetPhase = pi/2;

modes = {'fast','thyme'};
for iMode = 1:numel(modes)
  Opt.SimulationMode = modes{iMode};
  try
    saffron(Sys,Exp,Opt);
    ok(iMode) = false;
  catch err
    ok(iMode) = contains(err.message,'Exp.DetPhase cannot be used with predefined experiments');
  end
end
