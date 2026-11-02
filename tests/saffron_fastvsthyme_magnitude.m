function ok = test()

% Fast and thyme mode give the same signal amplitude Tr(S+*rho)/Tr(S+'*S+),
% for ideal and finite pulses. Thyme mode selects the echo pathway with
% phase cycles (single spin packet, no inhomogeneous broadening).

Sys.S = 1/2;
Sys.g = gfree;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.MolFrame = [0 0 0];
Exp.SampleFrame = [0 pi/4 0];

tau = 0.15;  % µs
T = 0.1;  % µs
dt = 0.02;  % µs

% Phase cycles: sel keeps Delta p = +1 (echo pathway after the first pulse),
% odd keeps odd and even keeps even Delta p
sel = [0 1; pi/2 1i; pi -1; 3*pi/2 -1i];
odd = [0 1; pi -1];
even = [0 1; pi 1];

tp = 0.001;  % µs, finite pulse length
p90 = struct('Flip',pi/2,'Phase',pi/2);
p180 = struct('Flip',pi,'Phase',pi/2);
f90 = p90; f90.tp = tp;
f180 = p180; f180.tp = 2*tp;

E2 = Exp; E2.Dim1 = {'d1,d2' dt}; E2.nPoints = 3;
E3 = Exp; E3.Dim1 = {'d2' dt}; E3.nPoints = 3;
EH = Exp; EH.Dim1 = {'d2' dt}; EH.Dim2 = {'d3' dt}; EH.nPoints = [2 2]; EH.DetPhase = pi;

E2.Sequence = {p90 tau p180 tau};
E3.Sequence = {p90 tau p90 T p90 tau};
EH.Sequence = {p90 tau p90 T p180 T p90 tau};
yIdeal = {fast(Sys,E2) fast(Sys,E3) fast(Sys,EH)};

E2.Sequence = {f90 tau f180 tau};
E3.Sequence = {f90 tau f90 T f90 tau};
EH.Sequence = {f90 tau f90 T f180 T f90 tau};
yFinite = {fast(Sys,E2) fast(Sys,E3) fast(Sys,EH)};

E2.PhaseCycle = {sel};
E3.PhaseCycle = {sel [] odd};
EH.PhaseCycle = {sel [] even odd};
yThyme = {thyme(Sys,E2) thyme(Sys,E3) thyme(Sys,EH)};

for k = 1:3
  ok(k) = areequal(yFinite{k},yThyme{k},1e-4,'rel');
  ok(3+k) = areequal(yIdeal{k},yThyme{k},5e-3,'rel');
end

% Without modulation, the stimulated and the four-pulse echo are half the
% primary echo
Sys.A = [0 0];
E2 = rmfield(E2,'PhaseCycle');
E3 = rmfield(E3,'PhaseCycle');
EH = rmfield(EH,'PhaseCycle');
E2.Sequence = {p90 tau p180 tau};
E3.Sequence = {p90 tau p90 T p90 tau};
EH.Sequence = {p90 tau p90 T p180 T p90 tau};
y2 = fast(Sys,E2);
ok(7) = areequal(fast(Sys,E3),y2/2,1e-10,'rel');
ok(8) = areequal(fast(Sys,EH),y2(1)/2*ones(4,1),1e-10,'rel');

end

function y = fast(Sys,Exp)
Exp.ExciteWidth = Inf;
Opt.SimulationMode = 'fast';
Opt.nOffsets = 1;
Opt.TimeDomain = true;  % exact time-domain signal instead of histogram
[~,y] = saffron(Sys,Exp,Opt);
y = y(:);
end

function y = thyme(Sys,Exp)
Exp.DetWindow = 0;
Opt.SimulationMode = 'thyme';
[~,y] = saffron(Sys,Exp,Opt);
y = y(:);
end
