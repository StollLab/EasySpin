function ok = test()

% Signal phases do not depend on the pulse timing or on the simulation
% frame: the microwave carrier is phase-continuous over the sequence

tol = 2;  % degrees
ok = [];

Sys.S = 1/2;
Sys.ZeemanFreq = 9.5;  % GHz

P90.Type = 'rectangular';
P90.tp = 0.010;  % µs
P90.Flip = pi/2;  % rad
P180 = P90;
P180.Flip = pi;

Exp.mwFreq = 9.5;  % GHz

% Two-pulse echo, transient and single-point detection: the echo phase is
% the same for sub-ns changes of tau and for different simulation frames
tau = 0.2 + [0 0.0003 0.0006];  % µs
SimFreq = [7 7.3];  % GHz
for detWindow = {[-0.002 0.002], 0}
  Exp.DetWindow = detWindow{1};
  phi = zeros(numel(SimFreq),numel(tau));
  for iFreq = 1:numel(SimFreq)
    Opt.SimFreq = SimFreq(iFreq);
    for iTau = 1:numel(tau)
      Exp.Sequence = {P90 tau(iTau) P180 tau(iTau)+P90.tp/2};
      [~,signal] = spidyan(Sys,Exp,Opt);
      [~,idx] = max(abs(signal));
      phi(iFreq,iTau) = angle(signal(idx));
    end
  end
  ok(end+1) = phasediff(phi,phi(1)) < tol; %#ok<AGROW>
end

% Single-point FID after a delay: same phase in different simulation frames
Exp.DetWindow = 0;
Exp.Sequence = {0.1003 P90 0.1171};
phi = zeros(1,3);
SimFreq = [7 7.3 7.5];
for iFreq = 1:numel(SimFreq)
  Opt.SimFreq = SimFreq(iFreq);
  [~,signal] = spidyan(Sys,Exp,Opt);
  phi(iFreq) = angle(signal);
end
ok(end+1) = phasediff(phi,phi(1)) < tol;

% Sz after two pi/2 pulses depends on their relative phase, which must not
% depend on the simulation frame
Exp.Sequence = {P90 0.1003 P90};
Exp.DetOperator = {'z1'};
Exp.DetFreq = 0;
Sz = zeros(1,2);
SimFreq = [7 7.3];
for iFreq = 1:numel(SimFreq)
  Opt.SimFreq = SimFreq(iFreq);
  [~,signal] = spidyan(Sys,Exp,Opt);
  Sz(iFreq) = real(signal);
end
ok(end+1) = abs(Sz(1)-Sz(2)) < 1e-2;

% Single points of a Hermitian detection operator (Sx) are not translated
% and stay real, while S+ is translated
Exp.Sequence = {P90 0.0537};
Exp.DetOperator = {'x1','+1'};
Exp.DetFreq = [9.5 9.5];
Sp = zeros(1,2);
for iFreq = 1:numel(SimFreq)
  Opt.SimFreq = SimFreq(iFreq);
  [~,signal] = spidyan(Sys,Exp,Opt);
  ok(end+1) = isreal(signal(1)) || abs(imag(signal(1))) < 1e-12; %#ok<AGROW>
  Sp(iFreq) = signal(2);
end
ok(end+1) = abs(Sp(2)-Sp(1)) < 0.02*abs(Sp(1));

% User-supplied IQ (lab frame, upconverted on the time axis of the pulse)
% acts like the built pulse, in different simulation frames
clear Opt
Opt.IntTimeStep = 4e-6;  % µs
Exp = rmfield(Exp,{'DetOperator','DetFreq'});
P45 = P90;
P45.Flip = pi/4;
P45.TimeStep = Opt.IntTimeStep;
[t,IQ] = pulse(P45);
OptMix.dt = Opt.IntTimeStep;
[t,IQ] = rfmixer(t,IQ,Exp.mwFreq,'IQshift',OptMix);
P45user.IQ = IQ;
P45user.t = t;
signalBuilt = zeros(1,2);
signalUser = zeros(1,2);
for iFreq = 1:numel(SimFreq)
  Opt.SimFreq = SimFreq(iFreq);
  Exp.Sequence = {P90 0.1003 P45 0.01};
  [~,signalBuilt(iFreq)] = spidyan(Sys,Exp,Opt);
  Exp.Sequence = {P90 0.1003 P45user 0.01};
  [~,signalUser(iFreq)] = spidyan(Sys,Exp,Opt);
end
ok(end+1) = phasediff(angle([signalBuilt signalUser]),angle(signalBuilt(1))) < tol && ...
  max(abs(abs(signalUser)-abs(signalBuilt))) < 0.02*abs(signalBuilt(1));

end

function d = phasediff(phi,phi0)
% largest phase difference in degrees, taking 2*pi periodicity into account
d = max(abs(rad2deg(angle(exp(1i*(phi(:)-phi0))))));
end
