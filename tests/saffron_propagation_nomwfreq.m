function ok = test() %#ok<FNDEF>

% Simulations without Exp.mwFreq: the microwave carrier is at 0, so pulse
% frequencies are absolute and the detection frame defaults to the lab
% frame. A pulse at 9500 MHz without Exp.mwFreq equals a pulse at offset 0
% with Exp.mwFreq = 9.5 GHz (same detection frame). In the full-space
% engine, finite pulses without a nonzero frequency and a resonator give
% errors without Exp.mwFreq. Pulses without a frequency use the subspace engine
% (on resonance), pulses with a frequency the full-space engine.

Sys.S = 1/2;
Sys.g = gfree*mhz2mt(9505)/mhz2mt(9500);  % 5 MHz above 9.5 GHz

Pulse.Type = 'rectangular';
Pulse.tp = 0.02;  % µs
Pulse.Flip = pi/2;

Exp.Field = mhz2mt(9500);  % mT
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = [0 0.05];  % µs
Exp.DetOperator = {'+1','z1'};
Exp.DetFreq = 9.5;  % GHz

% With carrier: offset 0
E = Exp;
E.mwFreq = 9.5;  % GHz
E.Sequence = {Pulse};
[~,yRef] = saffron(Sys,E);

% Without carrier: absolute frequency
E = Exp;
P = Pulse;
P.Frequency = 9500;  % MHz
E.Sequence = {P};
[~,y] = saffron(Sys,E);
ok(1) = areequal(y,yRef,1e-6,'abs');

% Default detection frame without Exp.mwFreq: lab frame (Sz unaffected)
E = rmfield(E,'DetFreq');
[~,yLab] = saffron(Sys,E);
ok(2) = areequal(yLab(:,2),yRef(:,2),1e-6,'abs');

% An empty Exp.mwFreq is the same as no Exp.mwFreq
E.mwFreq = [];
[~,yEmpty] = saffron(Sys,E);
ok(8) = areequal(yEmpty,yLab,1e-12,'abs');

% Errors
E = rmfield(Exp,'DetFreq');
E.Sequence = {Pulse};  % no frequency
ok(3) = throwsWith(@() saffron(Sys,E),'Pulse.Frequency');
P0 = Pulse;
P0.Frequency = 0;  % would be 0 MHz in the lab frame
E.Sequence = {P0};
ok(9) = throwsWith(@() saffron(Sys,E),'nonzero Pulse.Frequency');
E.Sequence = {P};
E.ResonatorFrequency = 9.5;  % GHz
E.ResonatorQL = 100;
ok(4) = throwsWith(@() saffron(Sys,E),'Exp.mwFreq');

% Pulses without a frequency stay on the subspace engine, which puts them
% on resonance with every orientation (ideal and finite pulses)
SysN.Nucs = '1H';
SysN.A = [2 3 9];  % MHz
ExpN.Field = 350;  % mT
ExpN.Sequence = '2pESEEM';
ExpN.tau = 0.1;  % µs
ExpN.dt = 0.02;  % µs
ExpN.nPoints = 16;
[~,~,info] = saffron(SysN,ExpN,struct('GridSize',7));
ok(5) = strcmp(info.Method,'subspace');
ExpN.tp = 0.016;  % µs, finite macro pulses
[~,~,info] = saffron(SysN,ExpN,struct('GridSize',7));
ok(6) = strcmp(info.Method,'subspace');

% A finite pulse with a frequency (absolute without Exp.mwFreq) needs the
% full-space engine
E = rmfield(Exp,{'DetOperator','DetFreq'});
E.DetWindow = 0;
E.Sequence = {P};
[~,~,info] = saffron(SysN,setfield(E,'Field',ExpN.Field),struct('GridSize',7)); %#ok<SFLD>
ok(7) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'absolute pulse frequencies');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
