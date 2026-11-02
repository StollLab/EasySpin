function ok = test()

% Signal phases do not depend on the pulse timing or on the simulation
% frame, and integrated echoes equal the integral of the transient

tol = 0.5;  % degrees

Sys.S = 1/2;
Sys.g = gfree;

Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);

P90.Type = 'rectangular';
P90.tp = 0.010;  % µs
P90.Flip = pi/2;  % rad
P180 = P90;
P180.Flip = pi;

Exp.DetWindow = [-0.002 0.002];  % µs

Opt.GridSize = 7;
Opt.SimulationMode = 'thyme';

% Two-pulse echo: the echo phase is the same for sub-ns changes of tau
% and for different simulation frames
tau = 0.2 + [0 0.0003];  % µs
SimFreq = [7 7.3];  % GHz
phi = zeros(numel(SimFreq),numel(tau));
for iFreq = 1:numel(SimFreq)
  Opt.SimFreq = SimFreq(iFreq);
  for iTau = 1:numel(tau)
    Exp.Sequence = {P90 tau(iTau) P180 tau(iTau)+P90.tp/2};
    [~,signal] = saffron(Sys,Exp,Opt);
    [~,idx] = max(abs(signal));
    phi(iFreq,iTau) = angle(signal(idx));
  end
end
ok(1) = max(abs(rad2deg(angle(exp(1i*(phi(:)-phi(1))))))) < tol;

% Integrated echo equals the integral of the down-converted transient
% (which includes a factor of 2 for complex signals), and does not depend
% on the detection time step
[t,transient] = saffron(Sys,Exp,Opt);
Exp.DetIntegrate = true;
integrated = saffron(Sys,Exp,Opt);
ok(2) = areequal(integrated,trapz(t,transient)/2,1e-6,'rel');
Exp.DetTimeStep = 2*(t(2)-t(1));
integrated2 = saffron(Sys,Exp,Opt);
ok(3) = areequal(integrated2,integrated,1e-2,'rel');

end
