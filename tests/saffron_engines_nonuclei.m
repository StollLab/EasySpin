function ok = test() %#ok<FNDEF>

% Spin systems without nuclei in the nuclear subspace
% - no nuclei: Opt.Method = 'auto' uses the manifold engine (each EPR
%   transition is a two-level system), which agrees with the full
%   engine for ideal and finite pulses
% - all nuclei in the core system (a strongly coupled nucleus moved there by
%   Opt.IdealPulses): no error, and an echo without modulation
% - the resonance offset enters the free evolution: FIDs after ideal and
%   finite pulses, and the unrefocused pathways of a two-pulse sequence
%   with a nucleus, agree between the engines
% - with ideal pulses, an FID decays with the inhomogeneous line
%   (Sys.lwAnisoG), and with Opt.IdealPulses with the excited part of it
%   (product of line and window)
% - errors: RF pulses without nuclei

Sys.g = [2.00 2.01 2.03];
Exp.mwFreq = 9.5;  % GHz
Exp.Field = 338;  % mT
Exp.MolFrame = [0.2 0.8 0.1];
Exp.CrystalSymmetry = 1;

% Echo pathway selection by phase cycles in the full engine: the
% first pulse in four steps of 90°, and for the stimulated echo also the
% third pulse in two steps of 180°, with the receiver phase following
ph = (0:3).'*pi/2;
sel2 = [ph 0*ph -ph];
[ph1,ph3] = ndgrid((0:3)*pi/2,[0 pi]);
sel3 = [ph1(:) 0*ph1(:) ph3(:) -ph1(:)-ph3(:)];
p90.Flip = pi/2;
p180.Flip = pi;
f90 = p90; f90.tp = 0.01;  % µs
f180 = p180; f180.tp = 0.02;  % µs

E2 = Exp; E2.Dim1 = {'d1,d2' 0.05}; E2.nPoints = 8; E2.DetPhase = -pi/2;
E3 = Exp; E3.Dim1 = {'d2' 0.05}; E3.nPoints = 8; E3.DetPhase = -pi/2;
Exps = {E2 E2 E3};
Exps{1}.Sequence = {p90 0.2 p180 0.2};
Exps{2}.Sequence = {f90 0.2 f180 0.2};
Exps{3}.Sequence = {f90 0.2 f90 0.1 f90 0.2};
PhaseCycles = {sel2 sel2 sel3};
tol = [1e-10 1e-3 1e-3];  % finite pulses: time step of the full engine

for k = 1:3
  E = Exps{k};
  E.Pathways = 'echo';
  [~,yMan,info] = saffron(Sys,E);
  E = Exps{k};
  E.PhaseCycle = PhaseCycles{k};
  E.DetWindow = 0;
  [~,yFull] = saffron(Sys,E,struct('Method','full'));
  ok(2*k-1) = strcmp(info.Method,'manifold'); %#ok<AGROW>
  ok(2*k) = areequal(yMan,yFull,tol(k),'rel'); %#ok<AGROW>
end

% Primary echo with ideal pulses: 1/2, no modulation
ok(end+1) = areequal(real(yManIdeal(Sys,Exps{1})),0.5*ones(1,8),1e-10,'abs');

% All nuclei in the core system
SysC = Sys;
SysC.Nucs = '63Cu';
SysC.A = [50 50 500];  % MHz, larger than the excitation window 1/tp = 50 MHz
E = Exps{1};
E.Sequence{1}.tp = 0.02;  % µs
E.Sequence{3}.tp = 0.04;  % µs
E.Pathways = 'echo';
[~,y,info] = saffron(SysC,E,struct('IdealPulses',true));
ok(end+1) = strcmp(info.Method,'manifold') && max(abs(y-y(1)))<1e-10*abs(y(1)) && abs(y(1))>0;

% FID 5 MHz off resonance, sampled at several delays after the pulse
SysF.g = gfree*mhz2mt(9505)/mhz2mt(9500);
ExpF.Field = mhz2mt(9500);  % mT
ExpF.mwFreq = 9.5;  % GHz
ExpF.MolFrame = [0 0 0];
ExpF.Dim1 = {'d1' 0.01};  % µs
ExpF.nPoints = 6;
for tp = [0 0.02]  % ideal and finite pulse
  p.Flip = pi/2;
  p.tp = tp;
  ExpF.Sequence = {p 0};
  [~,yMan] = saffron(SysF,ExpF,struct('Method','manifold'));
  E = ExpF;
  E.DetWindow = 0;
  [~,yFull] = saffron(SysF,E,struct('Method','full'));
  ok(end+1) = areequal(yMan,yFull,1e-3,'rel'); %#ok<AGROW>
end

% All pathways of a two-pulse sequence with a nucleus, off resonance
SysN = SysF;
SysN.Nucs = '1H';
SysN.A = [1 2 4];  % MHz
E = ExpF;
E.Sequence = {f90 0.1 f180 0.1};
E.Dim1 = {'d1,d2' 0.02};
E.nPoints = 8;
[~,yMan] = saffron(SysN,E,struct('Method','manifold'));
E.DetWindow = 0;
[~,yFull] = saffron(SysN,E,struct('Method','full'));
ok(end+1) = areequal(yMan,yFull,1e-3,'rel');

% FID on resonance with an inhomogeneous line of FWHM lw: Gaussian decay
% exp(-(pi*w*t)^2/(4*ln2)), with w = lw for ideal pulses, and w =
% lw*W/sqrt(lw^2+W^2) for an excitation window of FWHM W = 1/tp.
lw = 10;  % MHz
SysL.g = gfree;
SysL.lwAnisoG = lw;
E = ExpF;
E.Field = mhz2mt(9500);  % on resonance
E.Dim1 = {'d1' 0.02};
E.nPoints = 6;
t = (0:E.nPoints-1)*0.02;  % µs
p = struct('Flip',pi/2,'tp',0);
E.Sequence = {p 0};
[~,y] = saffron(SysL,E,struct('Method','manifold'));
ok(end+1) = areequal(abs(y)/abs(y(1)),exp(-(pi*lw*t).^2/(4*log(2))),1e-5,'abs');
p.tp = 0.05;  % µs, window W = 20 MHz
E.Sequence = {p 0};
[~,y] = saffron(SysL,E,struct('Method','manifold','IdealPulses',true));
W = 1/p.tp;
w = lw*W/sqrt(lw^2+W^2);
ok(end+1) = areequal(abs(y)/abs(y(1)),exp(-(pi*w*t).^2/(4*log(2))),1e-5,'abs');

% Errors
E = ExpF;
rf = struct('Channel','rf','tp',1,'Flip',pi);
E.Sequence = {p90 0.2 rf 0.1 p90 0.2};
E.rfFreq = 14;  % MHz
ok(end+1) = throwsWith(@() saffron(SysF,E),'no nuclei');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch ME
  ok = contains(ME.message,str);
end
end

function y = yManIdeal(Sys,Exp)
Exp.Pathways = 'echo';
[~,y] = saffron(Sys,Exp,struct('Method','manifold'));
end
