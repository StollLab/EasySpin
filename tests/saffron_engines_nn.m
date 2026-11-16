function ok = test() %#ok<FNDEF>

% Nuclear-nuclear couplings (Sys.nn) in both engines
% - single orientation: the subspace and full-space engines agree for 2p and
%   3p ESEEM with ideal and finite pulses, and nn changes the signal
% - Sys.nn = 0 gives the same result as no Sys.nn
% - nn makes the echo of a mirror-symmetric sequence complex at a single
%   orientation (also an isotropic J), so a partial-sphere powder grid needs
%   the field-reversed orientations; with them, it agrees with the full sphere
% - errors: product rule, nn-coupled nuclei outside Opt.Nuclei, and
%   strongly coupled nuclei with nn and Opt.IdealPulses

Sys.g = gfree;
Sys.Nucs = '1H,1H';
Sys.A = [0.5 1 2; -1 -0.5 1.5];  % MHz
Sys.AFrame = [0 0 0; 0.4 0.9 0.2];
Sys.nn = [0.3 0.4 -0.7];  % MHz
Sys.nnFrame = [0.2 0.6 0];

Sys0 = Sys;
Sys0.nn = [0 0 0];
SysNoNN = rmfield(Sys,{'nn','nnFrame'});

Exp.mwFreq = 9.8;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.MolFrame = [0.3 0.7 0.1];
Exp.CrystalSymmetry = 1;

dt = 0.05;  % µs
tau = 0.15;  % µs
T = 0.1;  % µs
p90.Flip = pi/2;
p180.Flip = pi;
f90 = p90; f90.tp = 0.002;  % µs
f180 = p180; f180.tp = 0.004;  % µs

% Echo pathway selection by phase cycles in the full-space engine: the
% first pulse in four steps of 90°, and for the stimulated echo also the
% third pulse in two steps of 180°, with the receiver phase following
ph = (0:3).'*pi/2;
sel2 = [ph 0*ph -ph];
[ph1,ph3] = ndgrid((0:3)*pi/2,[0 pi]);
sel3 = [ph1(:) 0*ph1(:) ph3(:) -ph1(:)-ph3(:)];

E2 = Exp; E2.Dim1 = {'d1,d2' dt}; E2.nPoints = 32; E2.DetPhase = -pi/2;
E3 = Exp; E3.Dim1 = {'d2' dt}; E3.nPoints = 32; E3.DetPhase = -pi/2;
Exps = {E2 E3 E2 E3};
Exps{1}.Sequence = {p90 tau p180 tau};
Exps{2}.Sequence = {p90 tau p90 T p90 tau};
Exps{3}.Sequence = {f90 tau f180 tau};
Exps{4}.Sequence = {f90 tau f90 T f90 tau};
PhaseCycles = {sel2 sel3 sel2 sel3};
% ideal pulses: high-field approximation of the subspace engine (A^2/nu_e);
% finite pulses: also the time step of the full-space engine
tol = [1e-5 1e-5 1e-4 1e-4];

for k = 1:4
  ySub = subspace(Sys,Exps{k});
  E = Exps{k};
  E.PhaseCycle = PhaseCycles{k};
  E.DetWindow = 0;
  [~,yFull] = saffron(Sys,E,struct('Method','fullspace'));
  ok(k) = areequal(ySub,yFull(:),tol(k),'rel'); %#ok<AGROW>
  % the couplings change the signal much more than the tolerance
  ok(4+k) = norm(ySub-subspace(Sys0,Exps{k}))/norm(ySub) > 1e-3; %#ok<AGROW>
end

% Sys.nn = 0 is the same as no Sys.nn
ok(end+1) = isequal(subspace(Sys0,Exps{3}),subspace(SysNoNN,Exps{3}));

% Complex echo at a single orientation: with nn (anisotropic, or isotropic
% J), but not without
SysJ = SysNoNN;
SysJ.nn = 0.5;  % MHz
yJ = subspace(SysJ,Exps{1});
yNN = subspace(Sys,Exps{1});
y0 = subspace(Sys0,Exps{1});
ok(end+1) = norm(imag(yJ))/norm(yJ) > 1e-5;
ok(end+1) = norm(imag(yNN))/norm(yNN) > 1e-5;
ok(end+1) = norm(imag(y0))/norm(y0) < 1e-12;

% Powders: partial-sphere grid (with field-reversed orientations) vs full
% sphere, for the conj shortcut (2p ESEEM with ideal pulses) and computed
% partners (finite pulses)
Powder = rmfield(Exp,{'MolFrame','CrystalSymmetry'});
EP{1} = Powder; EP{1}.Sequence = '2pESEEM'; EP{1}.dt = dt; EP{1}.nPoints = 16;
EP{2} = Powder; EP{2}.Sequence = {f90 tau f180 tau}; EP{2}.Dim1 = {'d1,d2' dt}; EP{2}.nPoints = 16;
Opt.Method = 'subspace';
Opt.GridSize = 7;
Opt.GridFrame = [0 0 0];
for S = {SysJ Sys}
  for k = 1:2
    Opt.GridSymmetry = 'Ci';
    [~,yPartial] = saffron(S{1},EP{k},Opt);
    Opt.GridSymmetry = 'C1';
    [~,yFull] = saffron(S{1},EP{k},Opt);
    ok(end+1) = areequal(yPartial,yFull,1e-10,'rel'); %#ok<AGROW>
  end
end

% Errors
E = Exps{1};
E.Pathways = 'echo';
ok(end+1) = throwsWith(@() saffron(Sys,E,struct('Method','subspace','ProductRule',true)),'ProductRule');
ok(end+1) = throwsWith(@() saffron(Sys,E,struct('Method','subspace','Nuclei',1)),'Opt.Nuclei');
% a strongly coupled nucleus cannot be moved to the core system with
% Opt.IdealPulses if it has an nn coupling
SysS = Sys;
SysS.A(1,:) = [80 90 100];  % MHz, larger than the excitation window 1/tp = 50 MHz
E = Exps{3};
E.Sequence{1}.tp = 0.02;  % µs
E.Sequence{3}.tp = 0.04;  % µs
E.Pathways = 'echo';
ok(end+1) = throwsWith(@() saffron(SysS,E,struct('Method','subspace','IdealPulses',true)),'Opt.IdealPulses');

end

function y = subspace(Sys,Exp)
Exp.Pathways = 'echo';
Opt.Method = 'subspace';
[~,y] = saffron(Sys,Exp,Opt);
y = y(:);
end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch ME
  ok = contains(ME.message,str);
end
end
