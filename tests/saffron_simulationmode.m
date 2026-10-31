function ok = test()

% Choice of Opt.SimulationMode: fast if possible, error if fast is requested
% but not possible

Sys.S = 1/2;
Sys.g = 2;

Exp.Field = 350;
Exp.mwFreq = 9.8;
Exp.Sequence = '2pESEEM';
Exp.tau = 0.1;
Exp.dt = 0.02;
Exp.nPoints = 5;
Exp.ExciteWidth = Inf;

Opt.GridSize = 7;

% Fast mode is used by default and when requested
[~,y0] = saffron(Sys,Exp,Opt);
Opt.SimulationMode = 'fast';
[~,y1] = saffron(Sys,Exp,Opt);
ok(1) = isequal(y0,y1);

% Exp.DetWindow requires thyme mode: used by default, error if fast is requested
Exp = rmfield(Exp,'ExciteWidth');
Exp.DetWindow = 0;
try
  saffron(Sys,Exp,Opt);
  ok(2) = false;
catch err
  ok(2) = contains(err.message,'SimulationMode');
end
Opt = rmfield(Opt,'SimulationMode');
[~,y2] = saffron(Sys,Exp,Opt);
Opt.SimulationMode = 'thyme';
[~,y3] = saffron(Sys,Exp,Opt);
ok(3) = isequal(y2,y3);

% Invalid simulation mode
Opt.SimulationMode = 'slow';
try
  saffron(Sys,Exp,Opt);
  ok(4) = false;
catch
  ok(4) = true;
end
