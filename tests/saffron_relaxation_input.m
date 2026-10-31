function ok = test()

% Handling of Sys.T1, Sys.T2 and Sys.T1T2 input

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

% Sys.T1T2 is rejected in both simulation modes
SysR = Sys;
SysR.T1 = 1;
SysR.T1T2 = [1 1];
ok(1) = throws(@() saffron(SysR,Exp,Opt));
Opt.SimulationMode = 'thyme';
ok(2) = throws(@() saffron(SysR,Exp,Opt));

% Zero relaxation times are ignored in thyme mode
SysR = Sys;
SysR.T1 = 0;
SysR.T2 = 0;
[~,y0] = saffron(Sys,Exp,Opt);
[~,y] = saffron(SysR,Exp,Opt);
ok(3) = areequal(y,y0,1e-10,'rel');

% Matrix relaxation times select thyme mode by default and error in fast mode
SysR = Sys;
SysR.T2 = [0 1; 0 0];
Opt = rmfield(Opt,'SimulationMode');
ok(4) = ~throws(@() saffron(SysR,Exp,Opt));
Opt.SimulationMode = 'fast';
ok(5) = throws(@() saffron(SysR,Exp,Opt));

end

function threw = throws(fcn)
try
  [~,~] = fcn();
  threw = false;
catch
  threw = true;
end
end
