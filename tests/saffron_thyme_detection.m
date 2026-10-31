function ok = test()

% Detection fields in thyme mode: Exp.DetFreq is used for down-conversion,
% spidyan-specific detection fields are rejected

Sys.S = 1/2;
Sys.g = 2;

Exp.Field = 350;
Exp.mwFreq = 9.8;

Pulse.Type = 'rectangular';
Pulse.tp = 0.016;  % µs
Pulse.Flip = pi/2;  % rad
Exp.Sequence = {Pulse 0.1};
Exp.DetWindow = [0 0.05];  % µs

Opt.GridSize = 7;
Opt.SimulationMode = 'thyme';

% Exp.mwFreq is the default detection frequency
[~,y0] = saffron(Sys,Exp,Opt);
Exp.DetFreq = Exp.mwFreq;
[~,y1] = saffron(Sys,Exp,Opt);
ok(1) = areequal(y1,y0,1e-10,'rel');

% A different detection frequency shifts the signal by the frequency difference
Exp.DetFreq = Exp.mwFreq + 0.01;  % GHz
[t,y2] = saffron(Sys,Exp,Opt);
ratio = y2./y0;
ok(2) = areequal(ratio/ratio(1),exp(-2i*pi*10*(t-t(1))),1e-8,'abs');

% Exp.DetOperator and Exp.DetSequence are specific to spidyan
ExpD = Exp;
ExpD.DetOperator = {'z1'};
ok(3) = throwsWith(@() saffron(Sys,ExpD,Opt),'DetOperator');
ExpD = rmfield(Exp,'DetWindow');
ExpD.DetSequence = 'last';
ok(4) = throwsWith(@() saffron(Sys,ExpD,Opt),'DetSequence');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
