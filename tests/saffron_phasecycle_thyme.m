function ok = test()

% Exp.PhaseCycle selects the thyme method, and is rejected with an
% explicit fast method

Sys.S = 1/2;
Sys.g = gfree;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

p90 = struct('tp',0.016,'Flip',pi/2,'Phase',pi/2);
p180 = struct('tp',0.032,'Flip',pi,'Phase',pi/2);

Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.MolFrame = [0 0 0];
Exp.SampleFrame = [0 pi/4 0];
Exp.Sequence = {p90 0.1 p180 0.1};
Exp.Dim1 = {'d1,d2' 0.02};
Exp.nPoints = 3;
Exp.ExciteWidth = Inf;

% without phase cycle, the sequence runs in fast mode
[~,y] = saffron(Sys,Exp);
ok(1) = numel(y)==3;

% with phase cycle, thyme mode is selected, which needs Exp.DetWindow
Exp.PhaseCycle = {[0 1; pi/2 1i; pi -1; 3*pi/2 -1i]};
try
  saffron(Sys,Exp);
  ok(2) = false;
catch err
  ok(2) = contains(err.message,'Exp.DetWindow');
end

% explicit fast method
Opt.SimulationMode = 'fast';
try
  saffron(Sys,Exp,Opt);
  ok(3) = false;
catch err
  ok(3) = contains(err.message,'Exp.PhaseCycle');
end
