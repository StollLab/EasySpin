function ok = test() %#ok<FNDEF>

% Full engine, RF pulses in a shifted simulation frame: they are
% propagated in the nuclear RF frame, with the rotating-wave approximation
% for the RF field. With Opt.SimFreq = 0, they are stepped in the lab
% frame. The two agree within the error of the rotating-wave approximation
% (Bloch-Siegert shift, about nu2/(4*nuRF) relative to the ENDOR effect).
% Where this error exceeds 1%, RF pulses are stepped in the shifted frame.
% The RF operator is the nuclear Zeeman interaction: with the hyperfine
% enhancement, the lab-frame RF pulses would include the electron Zeeman
% interaction and need much finer steps (see saffron_endor_rfenhancement).

p90.Flip = pi/2;
Opt.Method = 'full';
Opt.rfEnhancement = false;
OptLab = Opt;
OptLab.SimFreq = 0;
tolRWA = 5e-3;

% Mims ENDOR, 1H, single orientation, RF carrier sweep
Sys.Nucs = '1H';
Sys.A = [1 1 3];  % MHz
rf.Channel = 'rf';
rf.tp = 10;  % µs
rf.B2 = 0.5/(rf.tp*larmorfrq('1H',1));  % mT, pi pulse for a bare proton
Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0.3 0];
Exp.CrystalSymmetry = 1;
Exp.rfFreq = 13;  % MHz
Exp.Sequence = {p90 0.2 p90 0 rf 0.1 p90 0.2};
Exp.Pathways = '+00-';
Exp.DetPhase = -pi/2;
Exp.rfPhaseLocked = true;  % one RF phase, keeps the lab-frame runs short
Exp.nPoints = 41;
Exp.Dim1 = {'rfFreq',0.1};
[~,y] = saffron(Sys,Exp,Opt);
[~,yLab] = saffron(Sys,Exp,OptLab);
ok(1) = areequal(y,yLab,tolRWA*max(abs(yLab)),'abs') && max(abs(yLab-yLab(1)))>0.3*max(abs(yLab));

% 14N with a hyperfine coupling larger than twice the nuclear Zeeman
% frequency: the nucleus precesses in opposite senses in the two electron
% manifolds
SysN.Nucs = '14N';
SysN.A = [18 18 24];  % MHz
SysN.Q = 0.3;  % MHz
rfN = rf;
rfN.B2 = 0.5/(rfN.tp*larmorfrq('14N',1));  % mT
ExpN = Exp;
ExpN.Sequence{5} = rfN;
ExpN.rfFreq = 7;  % MHz
ExpN.nPoints = 61;
[~,yN] = saffron(SysN,ExpN,Opt);
[~,yNLab] = saffron(SysN,ExpN,OptLab);
ok(2) = areequal(yN,yNLab,tolRWA*max(abs(yNLab)),'abs') && max(abs(yNLab-yNLab(1)))>0.1*max(abs(yNLab));

% A rectangular RF pulse at a fixed frequency is a single step in the RF
% frame, exact for any length: a pi pulse equals two halves back to back,
% at the carrier (phase-continuous) or at an offset
ExpS = Exp;
ExpS.nPoints = 5;
ExpS.Dim1 = {'rfFreq',0.5};
for offset = [0 0.3]
  rfS = rf;
  rfS.Frequency = offset;
  ExpS.Sequence{5} = rfS;
  [~,yWhole] = saffron(Sys,ExpS,Opt);
  rfHalf = rfS;
  rfHalf.tp = rfS.tp/2;
  rfHalf2 = rfHalf;
  rfHalf2.Phase = 2*pi*offset*rfHalf.tp;  % the offset waveform restarts at each pulse
  ExpS2 = ExpS;
  ExpS2.Sequence = [ExpS.Sequence(1:4) {rfHalf 0 rfHalf2} ExpS.Sequence(6:end)];
  ExpS2.Pathways = '+000-';
  [~,yHalves] = saffron(Sys,ExpS2,Opt);
  ok(end+1) = areequal(yHalves,yWhole,1e-10,'rel'); %#ok<AGROW>
end

% Phase of the RF field: a pi/2 pulse with phase pi/2 creates nuclear
% coherence, detected as Ix and Iy. The direction depends on the sense of
% precession, so it checks the sign of the nuclear phase. This needs an RF
% carrier with the same phase in every shot.
SysP.Nucs = '1H';
SysP.A = [0 0 20];  % MHz
SysP.initState = sop([1/2 1/2],'ez');
ExpP.Field = 350;  % mT
ExpP.mwFreq = 9.8;  % GHz
ExpP.MolFrame = [0 0 0];
ExpP.CrystalSymmetry = 1;
ExpP.rfPhaseLocked = true;
ExpP.DetOperator = {'ex','ey'};
ExpP.rfFreq = larmorfrq('1H',ExpP.Field) + SysP.A(3)/2;  % MHz, one NMR transition
rfP.Channel = 'rf';
rfP.tp = 5;  % µs
rfP.Flip = pi/2;
rfP.Phase = pi/2;
ExpP.Sequence = {rfP 0 p90 0};
yP = saffron(SysP,ExpP,Opt);
yPLab = saffron(SysP,ExpP,OptLab);
rfP.Phase = 0;
ExpP.Sequence = {rfP 0 p90 0};
yP0 = saffron(SysP,ExpP,Opt);
ok(end+1) = areequal(yP,yPLab,tolRWA*max(abs(yPLab)),'abs') && ...
  ~areequal(yP,yP0,0.1*max(abs(yP0)),'abs');

% Relaxation during the RF pulse
SysR = Sys;
SysR.T1 = 20;  % µs
SysR.T2 = 3;  % µs
ExpR = Exp;
ExpR.nPoints = 21;
ExpR.Dim1 = {'rfFreq',0.2};
[~,yR] = saffron(SysR,ExpR,Opt);
[~,yRLab] = saffron(SysR,ExpR,OptLab);
ok(end+1) = areequal(yR,yRLab,tolRWA*max(abs(yRLab)),'abs') && ~areequal(yR,y(1:2:end),0.01*max(abs(y)),'abs');

% Detection during the RF pulse (whole sequence): the end of the transient
ExpD = rmfield(Exp,{'nPoints','Dim1','Pathways'});
ExpD.rfFreq = 13.4;  % MHz
ExpD.DetWindow = 'sequence';
[~,yD] = saffron(Sys,ExpD,Opt);
[~,yDLab] = saffron(Sys,ExpD,OptLab);
ok(end+1) = areequal(yD(end),yDLab(end),tolRWA*max(abs(yDLab)),'abs');

% Strong RF field at a low frequency: the rotating-wave approximation is
% inaccurate, so the RF pulse is stepped in the shifted microwave frame,
% without warning. This is exact: it agrees with the lab frame within the
% high-field approximation of the shifted frame, with a random RF phase
% (averaged until converged) and with a phase-locked RF carrier.
SysW.Nucs = '14N';
SysW.A = 0.5;  % MHz
rfW = rf;
rfW.tp = 1;  % µs
rfW.B2 = 100;  % mT
ExpW = Exp;
ExpW.Sequence{5} = rfW;
ExpW.rfFreq = 1.1;  % MHz
ExpW.nPoints = 3;
ExpW.Dim1 = {'rfFreq',0.3};
tolStep = 1e-5;
for locked = [false true]
  ExpW.rfPhaseLocked = locked;
  lastwarn('');
  [~,yW] = saffron(SysW,ExpW,Opt);
  noWarning = isempty(lastwarn);
  [~,yWLab] = saffron(SysW,ExpW,OptLab);
  ok(end+1) = noWarning && areequal(yW,yWLab,tolStep,'rel'); %#ok<AGROW>
end

% Stepped RF pulse with relaxation, and detected during the pulse
SysWR = SysW;
SysWR.T1 = 5;  % µs
SysWR.T2 = 0.5;  % µs
ExpWR = ExpW;
ExpWR.nPoints = 2;
[~,yWR] = saffron(SysWR,ExpWR,Opt);
[~,yWRLab] = saffron(SysWR,ExpWR,OptLab);
ExpWD = rmfield(ExpW,{'nPoints','Dim1','Pathways'});
ExpWD.DetWindow = 'sequence';
[~,yWD] = saffron(SysW,ExpWD,Opt);
[~,yWDLab] = saffron(SysW,ExpWD,OptLab);
ok(end+1) = areequal(yWR,yWRLab,tolStep,'rel') && ~areequal(yWR,yW(1:2),0.01,'rel');
% With detection, the stepped pulse uses the finer detection time grid,
% while the lab frame keeps at least 40 steps per RF period: they differ by
% the discretization error of the latter (about 1e-4 here)
ok(end+1) = areequal(yWD(end),yWDLab(end),1e-3*max(abs(yWDLab)),'abs');
