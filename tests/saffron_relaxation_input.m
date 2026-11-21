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

Opt.GridSize = 7;

% Sys.T1T2 is rejected by both engines
SysR = Sys;
SysR.T1 = 1;
SysR.T1T2 = [1 1];
ok(1) = throws(@() saffron(SysR,Exp,Opt));
Opt.Method = 'full';
ok(2) = throws(@() saffron(SysR,Exp,Opt));

% Zero relaxation times are ignored by the full engine
SysR = Sys;
SysR.T1 = 0;
SysR.T2 = 0;
[~,y0] = saffron(Sys,Exp,Opt);
[~,y] = saffron(SysR,Exp,Opt);
ok(3) = areequal(y,y0,1e-10,'rel');

% Matrix relaxation times select propagation by default and error with the manifold engine
SysR = Sys;
SysR.T2 = [0 1; 0 0];
Opt = rmfield(Opt,'Method');
ok(4) = ~throws(@() saffron(SysR,Exp,Opt));
Opt.Method = 'manifold';
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
