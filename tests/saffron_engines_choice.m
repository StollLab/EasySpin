function ok = test() %#ok<FNDEF>

% Automatic engine choice (Opt.Method = 'auto'): the subspace engine when it
% supports every feature, propagation otherwise; an explicit 'subspace'
% errors and names the unsupported feature.

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz

p90.Flip = pi/2;
chirp.Type = 'quartersin/linear';
chirp.tp = 0.1;  % µs
chirp.trise = 0.02;  % µs
chirp.Frequency = [-50 50];  % MHz
chirp.Flip = pi/2;

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0 0];
Exp.CrystalSymmetry = 1;

% ideal pulses, single point: pathways
E = Exp;
E.Sequence = {p90 0.1 p90 0.1};
[~,~,info] = saffron(Sys,E);
ok(1) = strcmp(info.Method,'subspace');

% shaped pulse: propagation
E = Exp;
E.Sequence = {chirp 0.2 p90 0.1};
[~,~,info] = saffron(Sys,E);
ok(2) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'shaped');

% transient detection: propagation
E = Exp;
E.Sequence = {p90 0.1 p90 0.1};
E.DetWindow = [-0.01 0.01];
[~,~,info] = saffron(Sys,E,struct('IntTimeStep',1e-4));
ok(3) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'DetWindow');

% explicit pathways with an unsupported feature: error naming it
E = Exp;
E.Sequence = {chirp 0.2 p90 0.1};
try
  saffron(Sys,E,struct('Method','subspace'));
  ok(4) = false;
catch err
  ok(4) = contains(err.message,'shaped');
end

% initial state: propagation
S = Sys;
S.initState = {kron([1 0; 0 0],eye(2)),'eigen'};
E = Exp;
E.Sequence = {p90 0.1 p90 0.1};
[~,~,info] = saffron(S,E);
ok(5) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'initState');
