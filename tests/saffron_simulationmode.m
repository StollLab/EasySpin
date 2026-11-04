function ok = test() %#ok<FNDEF>

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

% Incremented delays in an order the fast method does not support (dimension
% 2 before dimension 1): thyme is used by default, which needs a detection
% window, and the error gives the reason; error if fast is requested
p90 = struct('Flip',pi/2,'tp',0.016);
p180 = struct('Flip',pi,'tp',0.032);
ExpC.Field = 350;
ExpC.mwFreq = 9.8;
ExpC.Sequence = {p90 0.1 p90 0.1 p180 0.1 p90 0.1};
ExpC.Dim1 = {'d3' 0.02};
ExpC.Dim2 = {'d2' 0.02};
ExpC.nPoints = [2 2];
OptC.GridSize = 7;
ok(5) = throwsWith(@() saffron(Sys,ExpC,OptC),'DetWindow') && ...
  throwsWith(@() saffron(Sys,ExpC,OptC),'incremented delays');
OptC.SimulationMode = 'fast';
ok(6) = throwsWith(@() saffron(Sys,ExpC,OptC),'incremented delays');

% Adjacent pulses (no delay in between): error if fast is requested
ExpA = ExpC;
ExpA.Sequence = {p90 p180 0.1 p180 0.1};
ExpA.Dim1 = {'d1' 0.02};
ExpA = rmfield(ExpA,'Dim2');
ExpA.nPoints = 2;
ok(7) = throwsWith(@() saffron(Sys,ExpA,OptC),'not separated');

end

function ok = throwsWith(fcn,text)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,text);
end
end
