function ok = test() %#ok<FNDEF>

% Predefined experiments in thyme mode: initial delays default to zero as in
% fast mode, tau must be positive, and HYSCORE accepts one increment per
% dimension

Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.dt = 0.02;  % µs

Opt.GridSize = 7;
Opt.SimulationMode = 'thyme';

% 2pESEEM: tau defaults to 0
Exp2p = Exp;
Exp2p.Sequence = '2pESEEM';
Exp2p.nPoints = 3;
[~,yDefault] = saffron(Sys,Exp2p,Opt);
Exp2p.tau = 0;
[~,yZero] = saffron(Sys,Exp2p,Opt);
ok(1) = areequal(yDefault,yZero,1e-12,'rel');

% 3pESEEM: T defaults to 0
Exp3p = Exp;
Exp3p.Sequence = '3pESEEM';
Exp3p.tau = 0.1;  % µs
Exp3p.nPoints = 3;
[~,yDefault] = saffron(Sys,Exp3p,Opt);
Exp3p.T = 0;
[~,yZero] = saffron(Sys,Exp3p,Opt);
ok(2) = areequal(yDefault,yZero,1e-12,'rel');

% HYSCORE: t1 and t2 default to 0
ExpHY = Exp;
ExpHY.Sequence = 'HYSCORE';
ExpHY.tau = 0.1;  % µs
ExpHY.nPoints = [2 2];
[~,yDefault] = saffron(Sys,ExpHY,Opt);
ExpHY.t1 = 0;
ExpHY.t2 = 0;
[~,yZero] = saffron(Sys,ExpHY,Opt);
ok(3) = areequal(yDefault,yZero,1e-12,'rel');

% tau must be larger than zero for 3pESEEM
Exp3p.tau = 0;
try
  saffron(Sys,Exp3p,Opt);
  ok(4) = false;
catch err
  ok(4) = contains(err.message,'Exp.tau must be larger than 0');
end

% HYSCORE with one increment per dimension, same as user-defined sequence
ExpHY.t1 = 0.05;  % µs
ExpHY.t2 = 0.07;  % µs
ExpHY.dt = [0.02 0.03];  % µs
ExpHY.nPoints = [3 2];
[x,yPredefined] = saffron(Sys,ExpHY,Opt);
P90 = struct('tp',0.016,'Flip',pi/2,'Phase',pi/2);
P180 = struct('tp',0.032,'Flip',pi,'Phase',pi/2);
ExpCustom = rmfield(ExpHY,{'tau','t1','t2','dt'});
ExpCustom.Sequence = {P90 ExpHY.tau P90 ExpHY.t1 P180 ExpHY.t2 P90 ExpHY.tau};
ExpCustom.Dim1 = {'d2' ExpHY.dt(1)};
ExpCustom.Dim2 = {'d3' ExpHY.dt(2)};
ExpCustom.DetWindow = 0;
ExpCustom.DetPhase = pi;
sel = [0 1; pi/2 1i; pi -1; 3*pi/2 -1i];
ExpCustom.PhaseCycle = {sel [] [0 1; pi 1] [0 1; pi -1]};
[~,yCustom] = saffron(Sys,ExpCustom,Opt);
ok(5) = areequal(yCustom,yPredefined,1e-10,'rel');
ok(6) = areequal(x{1},[0.05 0.07 0.09],1e-12,'abs') && areequal(x{2},[0.07 0.10],1e-12,'abs');

% 4pESEEM: T defaults to 0
Exp4p = Exp;
Exp4p.Sequence = '4pESEEM';
Exp4p.tau = 0.1;  % µs
Exp4p.nPoints = 3;
[~,yDefault] = saffron(Sys,Exp4p,Opt);
Exp4p.T = 0;
[~,yZero] = saffron(Sys,Exp4p,Opt);
ok(7) = areequal(yDefault,yZero,1e-12,'rel');