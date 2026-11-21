function ok = test() %#ok<FNDEF>

% Full engine, random RF phase (Exp.rfPhaseLocked = false, default): the
% RF carrier has a random phase in each shot, shared by all RF pulses of
% the shot. The signal is the average over this phase. It equals an
% explicit phase cycle over the shared RF phase (Exp.PhaseCycle with
% enough steps, phase-locked carrier), for one and two RF pulses (with and
% without a microwave pulse between them), for I = 1/2 and I = 1, with
% relaxation, and for stepped RF pulses. A phase-locked carrier gives a
% different signal in a single orientation.

p90.Flip = pi/2;
Opt.Method = 'full';
tol = 1e-10;
N = 12;  % steps of the explicit phase cycle, more than needed
phi = 2*pi*(0:N-1).'/N;

% Mims ENDOR, 1H, single orientation
Sys.Nucs = '1H';
Sys.A = [1 1 3];  % MHz
rf.Channel = 'rf';
rf.tp = 5;  % µs
rf.B1 = 0.5/(rf.tp*larmorfrq('1H',1));  % mT, pi pulse for a bare proton
Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0.3 0];
Exp.CrystalSymmetry = 1;
Exp.rfFreq = 13;  % MHz
Exp.Sequence = {p90 0.2 p90 0 rf 0.1 p90 0.2};
Exp.Pathways = '+00-';
Exp.DetPhase = -pi/2;
Exp.nPoints = 11;
Exp.Dim1 = {'rfFreq',0.2};
[ok(1),ok(2)] = compareCycle(Sys,Exp,Opt,[zeros(N,2) phi zeros(N,2)],tol);

% Two RF pulses with an incremented delay between them: they share the
% random phase
rf2 = rf;
rf2.B1 = rf.B1/2;
Exp2 = Exp;
Exp2.rfFreq = 13.4;  % MHz
Exp2.Sequence = {p90 0.2 p90 0 rf2 0.3 rf2 0.1 p90 0.2};
Exp2.Pathways = '+000-';
Exp2.nPoints = 4;
Exp2.Dim1 = {'d3',0.1};
[ok(3),ok(4)] = compareCycle(Sys,Exp2,Opt,[zeros(N,2) phi phi zeros(N,2)],tol);

% 14N (I = 1): harmonics of the RF phase up to order 4
SysN.Nucs = '14N';
SysN.A = [18 18 24];  % MHz
SysN.Q = [-0.3 -0.2 0.5];  % MHz
rfN = rf;
rfN.tp = 10;  % µs
rfN.B1 = 0.5/(rfN.tp*larmorfrq('14N',1));  % mT
ExpN = Exp;
ExpN.Sequence{5} = rfN;
ExpN.rfFreq = 7;  % MHz
[ok(5),ok(6)] = compareCycle(SysN,ExpN,Opt,[zeros(N,2) phi zeros(N,2)],tol);

% Relaxation during the sequence
SysR = Sys;
SysR.T1 = 20;  % µs
SysR.T2 = 3;  % µs
ExpR = Exp;
ExpR.nPoints = 3;
[ok(7),ok(8)] = compareCycle(SysR,ExpR,Opt,[zeros(N,2) phi zeros(N,2)],tol);

% NMR transitions at low frequencies: the RF pulse is stepped without the
% rotating-wave approximation, and the counter-rotating part of the field
% adds higher harmonics of the RF phase. The number of phases is doubled
% until the average converges (to 1e-4 relative to the largest signal).
SysN2 = SysN;
SysN2.A = [5 5 7];  % MHz
ExpN2 = ExpN;
ExpN2.Sequence{5}.tp = 5;  % µs
ExpN2.Sequence{5}.B1 = 2*rfN.B1;
ExpN2.rfFreq = 2.5;  % MHz
M = 24;
[ok(9),ok(10)] = compareCycle(SysN2,ExpN2,Opt,[zeros(M,2) 2*pi*(0:M-1).'/M zeros(M,2)],1e-4);

% Two RF pulses with a microwave pi pulse between them: the RF pulses no
% longer act as one, and the signal has harmonics of the RF phase up to
% twice the order of a single RF pulse
p180.Flip = pi;
Exp3 = Exp2;
Exp3.Sequence = {p90 0.2 p90 0 rf2 0.1 p180 0.1 rf2 0.1 p90 0.2};
Exp3.Pathways = '+0000-';
Exp3.Dim1 = {'d4',0.1};
[ok(11),ok(12)] = compareCycle(Sys,Exp3,Opt,[zeros(N,2) phi zeros(N,1) phi zeros(N,2)],tol);

% Exp.rfPhaseLocked must be true or false
try
  saffron(Sys,setfield(Exp,'rfPhaseLocked',2),Opt); %#ok<SFLD>
  ok(13) = false;
catch err
  ok(13) = contains(err.message,'Exp.rfPhaseLocked must be true or false');
end

end

% Random RF phase against an explicit phase cycle over the RF phase with a
% phase-locked carrier, and against the phase-locked carrier without cycle
function [okCycle,okLocked] = compareCycle(Sys,Exp,Opt,PhaseCycle,tol)
[~,yRandom] = saffron(Sys,Exp,Opt);
Exp.rfPhaseLocked = true;
[~,yLocked] = saffron(Sys,Exp,Opt);
Exp.PhaseCycle = PhaseCycle;
[~,yCycle] = saffron(Sys,Exp,Opt);
okCycle = areequal(yRandom,yCycle,tol,'rel');
okLocked = ~areequal(yRandom,yLocked,1e-3,'rel');
end
