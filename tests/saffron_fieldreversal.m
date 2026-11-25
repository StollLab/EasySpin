function ok = test() %#ok<FNDEF>

% Manifold engine, powder: with a nuclear quadrupole interaction, the
% signals of the field directions n and -n differ (complex-conjugated
% eigenvectors). A partial-sphere grid (Ci) must give the same signal as the
% full sphere (C1), for the shortcuts found by fieldReversalShortcut (2pESEEM:
% real part, HYSCORE: conjugate transpose) and for orientations computed a
% second time with the field reversed (finite pulses, sequences that are
% not mirror symmetric, point-by-point dimensions that break the symmetry
% of 2pESEEM at other points: a delay increment outside the binning schemes,
% and non-linear delay values).

Sys.Nucs = '14N';
Sys.A = [1.6 1.9 2.6];  % MHz
Sys.AFrame = [0.3 0.5 0.2];
Sys.Q = [0.9 0.6];  % MHz

Base.Field = 324.9;  % mT

Opt.Method = 'manifold';
Opt.GridSize = 7;

p90.Flip = pi/2;
p180.Flip = pi;
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
% 2pESEEM as user-defined experiment, with a point-by-point first dimension
E2p = Base;
E2p.Sequence = {p90 0.12 p180 0.12};
E2p.Pathways = '+-';
E2p.DetPhase = -pi/2;
E2p.nPoints = [4 16];
E2p.Dim2 = {'d1,d2' dt};
E = E2p; E.Dim1 = {'d1' 0.03};  % both dimensions increment the first interval
Exps{5} = E;
E = E2p; E.Dim1 = {'d1' [0 0.03 0.1 0.25]};  % non-linear values
Exps{6} = E;

for k = 1:numel(Exps)
  Opt.GridSymmetry = 'Ci';
  Opt.GridFrame = [0 0 0];
  [~,yPartial] = saffron(Sys,Exps{k},Opt);
  Opt.GridSymmetry = 'C1';
  [~,yFull] = saffron(Sys,Exps{k},Opt);
  ok(k) = areequal(yPartial,yFull,1e-10,'rel');
end

% without the first dimension, the shortcut applies to this 2pESEEM, so
% the cases above test that the point-by-point dimensions disable it
E = E2p; E.nPoints = 16; E.Dim1 = E.Dim2; E = rmfield(E,'Dim2');
Opt.GridSymmetry = 'Ci';
[~,~,info] = saffron(Sys,E,Opt);
ok(end+1) = info.reversedField=="conj";

end
