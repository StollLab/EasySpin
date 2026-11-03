function ok = test() %#ok<FNDEF>

% Fast mode, powder: with a nuclear quadrupole interaction, the signals of
% the field directions n and -n differ (complex-conjugated eigenvectors).
% A partial-sphere grid (Ci) must give the same signal as the full sphere
% (C1), for the predefined shortcuts (2pESEEM: real part, HYSCORE:
% conjugate transpose) and for orientations computed a second time with
% the field reversed (finite pulses, sequences that are not mirror
% symmetric).

Sys.Nucs = '14N';
Sys.A = [1.6 1.9 2.6];  % MHz
Sys.AFrame = [0.3 0.5 0.2];
Sys.Q = [0.9 0.6];  % MHz

Base.Field = 324.9;  % mT

Opt.SimulationMode = 'fast';
Opt.TimeDomain = true;
Opt.GridSize = 7;
Opt.nOffsets = 5;

p90 = struct('Flip',pi/2,'Phase',pi/2);
p180 = struct('Flip',pi,'Phase',pi/2);
f90 = p90; f90.tp = 0.016;  % µs
f180 = p180; f180.tp = 0.032;  % µs
dt = 0.05;  % µs

E = Base; E.Sequence = '2pESEEM'; E.dt = dt; E.nPoints = 16;
Exps{1} = E;
E = Base; E.Sequence = 'HYSCORE'; E.tau = 0.15; E.dt = dt; E.nPoints = [8 8];
Exps{2} = E;
E = Base; E.Sequence = {f90 0.1 f180 0.1}; E.Dim1 = {'d1,d2' dt}; E.nPoints = 16;
Exps{3} = E;
E = Base; E.Sequence = {p90 0.15 p90 0.1 p180 0.1 p90 0.15}; E.Dim1 = {'d3' dt}; E.nPoints = 16;
Exps{4} = E;

for k = 1:numel(Exps)
  Opt.GridSymmetry = 'Ci';
  Opt.GridFrame = [0 0 0];
  [~,yPartial] = saffron(Sys,Exps{k},Opt);
  Opt.GridSymmetry = 'C1';
  [~,yFull] = saffron(Sys,Exps{k},Opt);
  ok(k) = areequal(yPartial,yFull,1e-10,'rel');
end

end
