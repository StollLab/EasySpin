function ok = test() %#ok<FNDEF>

% Subspace engine: average over the inhomogeneous EPR line (Sys.lwAnisoG)
% - ideal 90-180 echo of an isolated electron spin on resonance, detected
%   around and far from the echo: Gaussian echo shape
%   exp(-(pi*lw*t)^2/(4*ln2)) at net coherence time t. Far from the echo,
%   the signal stays zero (no aliased echoes from the offset grid).
% - convergence with Opt.OffsetTolerance (finite pulses): errors of the
%   order of the tolerance; info.nOffsets grows as it decreases
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
Exp.Sequence = {p90 tau p180 0};
Exp.nPoints = 64;
echo = @(t) exp(-(pi*lw*t).^2/(4*log(2)))/2;

% around the echo
dt = 0.1/63;  % µs
d2 = tau - 0.05 + (0:63)*dt;
Exp.Sequence{4} = d2(1);
Exp.Dim1 = {'d2' dt};
[~,y] = saffron(Sys,Exp,Opt);
ok(1) = areequal(abs(y),echo(d2-tau),1e-5,'abs');

% far from the echo: net coherence times 1 to 3 µs
dt = 2/63;  % µs
d2 = tau + 1 + (0:63)*dt;
Exp.Sequence{4} = d2(1);
Exp.Dim1 = {'d2' dt};
[~,y] = saffron(Sys,Exp,Opt);
ok(2) = max(abs(y))<1e-5;

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

% errors
ok(7) = throwsWith(@() saffron(SysN,Exp2,struct('nOffsets',31)),'OffsetTolerance');
ok(8) = throwsWith(@() saffron(SysN,Exp2,struct('OffsetTolerance',0)),'OffsetTolerance');
ok(9) = throwsWith(@() saffron(SysN,Exp2,struct('OffsetTolerance',1)),'OffsetTolerance');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch ME
  ok = contains(ME.message,str);
end
end
