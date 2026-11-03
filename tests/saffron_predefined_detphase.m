function ok = test() %#ok<FNDEF>

% Predefined experiments have a fixed detection phase: Exp.DetPhase is
% rejected in both simulation modes, also with a value of 0

Sys.Nucs = '1H';
Sys.A = 5;
Exp.Field = 340;
Exp.Sequence = '2pESEEM';
Exp.dt = 0.01;
Exp.DetPhase = 0;

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
