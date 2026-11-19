function ok = test() %#ok<FNDEF>

% Event-specific relaxation (Opt.Relaxation, one entry per element of
% Exp.Sequence) works with every detection window: the detection events
% appended for Exp.DetWindow take the setting of the last delay. The
% signal at the end of the sequence is the same for DetWindow = 0,
% [0 d] and 'sequence'.

Sys.S = 1/2;
Sys.T1 = 0.5;  % µs
Sys.T2 = 0.2;  % µs

Pulse.Type = 'rectangular';
Pulse.tp = 0.02;  % µs
Pulse.Flip = pi;

Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.mwFreq = 9.5;  % GHz
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetOperator = {'z1'};

Opt.Relaxation = [0 1];  % relaxation only during the delay
Opt.SimFreq = 9;  % GHz

% Single point at the end of the sequence
E = Exp;
E.Sequence = {Pulse 0.2};
E.DetWindow = 0;
y0 = saffron(Sys,E,Opt);

% Transient after the last delay: the window continues the relaxation of
% the last delay, so it ends at the same time with the same state
E = Exp;
E.Sequence = {Pulse 0.1};
E.DetWindow = [0 0.1];  % µs
yWin = saffron(Sys,E,Opt);

% Whole sequence
E = Exp;
E.Sequence = {Pulse 0.2};
E.DetWindow = 'sequence';
ySeq = saffron(Sys,E,Opt);

% <Sz> (normalized by Tr(Sz'*Sz)) relaxes from about +1 toward -1 during
% the delay
ok(1) = areequal(real(yWin(end)),real(y0),1e-3,'abs');
ok(2) = areequal(real(ySeq(end)),real(y0),1e-3,'abs');
ok(3) = real(yWin(end)) < real(yWin(1)) - 0.1;  % relaxation is on in the window

% Sequence ending with a pulse: the detection window after it takes the
% setting of the last delay, not of the pulse
E = Exp;
E.Sequence = {Pulse 0.1 Pulse};
E.DetWindow = [0 0.1];  % µs
OptP = Opt;
OptP.Relaxation = [0 1 0];
yPulseEnd = saffron(Sys,E,OptP);
E.Sequence = {Pulse 0.1 Pulse 0.1};
E.DetWindow = 0;
OptP.Relaxation = [0 1 0 1];
yRef = saffron(Sys,E,OptP);
ok(4) = areequal(real(yPulseEnd(end)),real(yRef),1e-3,'abs');

% A wrong number of entries is an error
E = Exp;
E.Sequence = {Pulse 0.2};
ok(5) = throwsWith(@() saffron(Sys,E,setfield(Opt,'Relaxation',[0 1 1])),'Opt.Relaxation');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
