function ok = test() %#ok<FNDEF>

% Subspace engine: average over the inhomogeneous EPR line (Sys.lwAnisoG)
% - ideal pulses, 90-180 echo of an isolated electron spin on resonance:
%   exact Gaussian echo shape exp(-(pi*lw*t)^2/(4*ln2)) at net coherence
%   time t, without offset grid
% - finite pulses, far from the echo: the signal stays zero (no aliased
%   echoes from the offset grid)
% - convergence with Opt.OffsetTolerance (finite pulses): errors of the
%   order of the tolerance; info.nOffsets grows as it decreases
% - ideal pulses, unrefocused pathways: the product rule gives the same
%   signal, and an ENDOR signal decays by exp(-(pi*lw*t)^2/(4*ln2))
% - errors: Opt.nOffsets, invalid Opt.OffsetTolerance

lw = 60;  % MHz
Sys.g = gfree;
Sys.lwAnisoG = lw;
Exp.mwFreq = 9.5;  % GHz
Exp.Field = mhz2mt(9500);  % on resonance
Exp.Pathways = '+-';
Opt.Method = 'subspace';
p90.Flip = pi/2;
p180.Flip = pi;
tau = 0.5;  % µs
Exp.nPoints = 64;
lineDecay = @(t,lw) exp(-(pi*lw*t).^2/(4*log(2)));

% ideal pulses, around the echo
dt = 0.1/63;  % µs
d2 = tau - 0.05 + (0:63)*dt;
Exp.Sequence = {p90 tau p180 d2(1)};
Exp.Dim1 = {'d2' dt};
[~,y,info] = saffron(Sys,Exp,Opt);
ok(1) = areequal(abs(y),lineDecay(d2-tau,lw)/2,1e-10,'abs') && info.nOffsets==1;

% finite pulses, far from the echo: net coherence times 1 to 3 µs
f90 = p90; f90.tp = 0.002;  % µs
f180 = p180; f180.tp = 0.004;  % µs
dt = 2/63;  % µs
d2 = tau + 1 + (0:63)*dt;
Exp.Sequence = {f90 tau f180 d2(1)};
Exp.Dim1 = {'d2' dt};
[~,y,info] = saffron(Sys,Exp,Opt);
ok(2) = max(abs(y))<1e-5 && info.nOffsets>1;

% convergence: finite pulses, nucleus, anisotropic line
SysN.g = gfree;
SysN.Nucs = '1H';
SysN.A = [2 6];  % MHz
SysN.lwAnisoG = [15 20 30];  % MHz
Exp2.mwFreq = 9.5;  % GHz
Exp2.Field = mhz2mt(9500);
Exp2.Sequence = '2pESEEM';
Exp2.tp = [0.016 0.032];  % µs
Exp2.dt = 0.02;  % µs
Exp2.nPoints = 128;
Opt2.GridSize = 10;
Opt2.OffsetTolerance = 1e-3;
[~,y3,info3] = saffron(SysN,Exp2,Opt2);
Opt2.OffsetTolerance = 1e-6;
[~,y6,info6] = saffron(SysN,Exp2,Opt2);
Opt2.OffsetTolerance = 1e-10;
[~,y10,info10] = saffron(SysN,Exp2,Opt2);
% errors of the order of the tolerance (truncation and aliasing, each up to tol)
ok(3) = areequal(y3,y10,5e-3,'rel');
ok(4) = areequal(y6,y10,5e-6,'rel');
ok(5) = info3.nOffsets<info6.nOffsets && info6.nOffsets<info10.nOffsets;

% without inhomogeneous line: a single offset
SysN.lwAnisoG = 0;
[~,~,info] = saffron(SysN,Exp2,Opt2);
ok(6) = info.nOffsets==1;

% ideal pulses, product rule: same signal, including the pathways that do
% not refocus (d1 incremented)
SysP.g = gfree;
SysP.Nucs = '1H,1H';
SysP.A = [1 2 3; 2 3 4];  % MHz
SysP.lwAnisoG = [15 20 30];  % MHz
Exp3.mwFreq = 9.5;  % GHz
Exp3.Field = mhz2mt(9500);
Exp3.Sequence = {p90 0.1 p180 0.1};
Exp3.Dim1 = {'d1' 0.02};
Exp3.nPoints = 4;
Opt3.Method = 'subspace';
Opt3.GridSize = 10;
Opt3.ProductRule = false;
[~,yNoProduct] = saffron(SysP,Exp3,Opt3);
Opt3.ProductRule = true;
[~,yProduct] = saffron(SysP,Exp3,Opt3);
ok(7) = areequal(yProduct,yNoProduct,1e-10,'rel');

% ideal pulses, ENDOR: the pathway '-00-' has net coherence time -2*tau
% and decays with the line (single orientation, isotropic line)
SysE.Nucs = '1H';
SysE.A = [1 1 3];  % MHz
rf = struct('Channel','rf','tp',5,'Flip',pi);
tau = 0.2;  % µs
ExpE.Field = 350;  % mT
ExpE.MolFrame = [0 0.3 0];
ExpE.CrystalSymmetry = 1;
ExpE.rfFreq = 12;  % MHz
ExpE.Sequence = {p90 tau p90 0 rf 0 p90 tau};
ExpE.Dim1 = {'rfFreq' 0.1};
ExpE.nPoints = 41;
ExpE.Pathways = '-00-';
[~,yNarrow] = saffron(SysE,ExpE,Opt);
SysE.lwAnisoG = 2;  % MHz
[~,yLine] = saffron(SysE,ExpE,Opt);
ok(8) = areequal(yLine,yNarrow*lineDecay(2*tau,2),1e-10,'rel');

% errors
ok(9) = throwsWith(@() saffron(SysN,Exp2,struct('nOffsets',31)),'OffsetTolerance');
ok(10) = throwsWith(@() saffron(SysN,Exp2,struct('OffsetTolerance',0)),'OffsetTolerance');
ok(11) = throwsWith(@() saffron(SysN,Exp2,struct('OffsetTolerance',1)),'OffsetTolerance');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch ME
  ok = contains(ME.message,str);
end
end
