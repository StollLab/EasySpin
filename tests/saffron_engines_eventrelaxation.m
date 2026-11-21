function ok = test() %#ok<FNDEF>

% Relaxation per element of Exp.Sequence (Opt.Relaxation) in the manifold
% engine: only the delays with relaxation on contribute their decay
% - stimulated echo: T2 decay over the two tau delays, T1 decay over T,
%   each only if switched on; entries for pulses have no effect
% - incremented delay T, evaluated per point (16 points) and binned (128
%   points): the T1 decay along T only if relaxation is on during T
% - two delays in one interval, only one of them with relaxation
% - Opt.Method = 'auto' uses the manifold engine, and agrees with the
%   full engine (T2 only, since the engines treat T1 differently)
% - a wrong number of entries is an error

Sys.Nucs = '1H';
Sys.A = 0;  % no modulation
Sys.T1 = 2;  % µs
Sys.T2 = 0.5;  % µs

p90.Flip = pi/2;
p180.Flip = pi;
tau = 0.15;  % µs
T = 0.4;  % µs

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0 0];
Exp.CrystalSymmetry = 1;
Exp.Sequence = {p90 tau p90 T p90 tau};
Exp.Pathways = '+0-';
Opt.Method = 'manifold';

% single point
y = @(relax) saffron(Sys,Exp,setfield(Opt,'Relaxation',relax));
y0 = y(false);
ok(1) = areequal(y(true),y0*exp(-2*tau/Sys.T2-T/Sys.T1),1e-12,'rel');
ok(2) = areequal(y([0 1 0 0 0 1]),y0*exp(-2*tau/Sys.T2),1e-12,'rel');
ok(3) = areequal(y([1 0 1 1 1 0]),y0*exp(-T/Sys.T1),1e-12,'rel');

% incremented T, per point and binned
dT = 0.02;  % µs
for nPts = [16 128]
  E = Exp;
  E.Dim1 = {'d2' dT};
  E.nPoints = nPts;
  Tk = T + (0:nPts-1)*dT;
  O = Opt;
  O.Relaxation = false;
  [~,yOff] = saffron(Sys,E,O);
  O.Relaxation = [0 1 0 0 0 1];
  [~,yTau] = saffron(Sys,E,O);
  O.Relaxation = [0 0 0 1 0 0];
  [~,yT] = saffron(Sys,E,O);
  ok(end+1) = areequal(yTau,yOff*exp(-2*tau/Sys.T2),1e-10,'rel'); %#ok<AGROW>
  ok(end+1) = areequal(yT,yOff.*exp(-Tk/Sys.T1),1e-10,'rel'); %#ok<AGROW>
end

% two delays in one interval: relaxation only during the first one
E = Exp;
E.Sequence = {p90 0.1 0.05 p180 0.15};
E.Pathways = '+-';
O = Opt;
O.Relaxation = false;
yOff = saffron(Sys,E,O);
O.Relaxation = [0 1 0 0 1];
ok(end+1) = areequal(saffron(Sys,E,O),yOff*exp(-0.25/Sys.T2),1e-12,'rel');

% automatic engine selection, and agreement with the full engine
SysT2 = rmfield(Sys,'T1');
E = Exp;
E.Sequence = {p90 0.1 p180 0.2};
E.Pathways = '+-';
E.DetPhase = -pi/2;
relax = [0 1 0 0];  % relaxation only during the first delay
[~,yMan,info] = saffron(SysT2,E,struct('Relaxation',relax));
[~,yFull] = saffron(SysT2,E,struct('Relaxation',relax,'Method','full'));
ok(end+1) = strcmp(info.Method,'manifold');
ok(end+1) = areequal(real(yFull),real(yMan),1e-3,'rel');
ok(end+1) = abs(yMan/saffron(SysT2,E,struct('Relaxation',false))-exp(-0.1/Sys.T2))<1e-12;

% wrong number of entries
try
  saffron(Sys,Exp,setfield(Opt,'Relaxation',[0 1 1]));
  ok(end+1) = false;
catch ME
  ok(end+1) = contains(ME.message,'Opt.Relaxation');
end
