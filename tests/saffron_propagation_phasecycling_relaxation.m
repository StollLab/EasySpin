function ok = test() %#ok<FNDEF>

orig_state = warning;
warning('off','all');

% Phase cycling with relaxation: the cycled signal must equal the weighted
% sum of separate experiments with the pulse phases of the cycle steps,
% including the recovery toward equilibrium after the cycled pulses
% (single orientation, detection after the last pulse, and over the whole
% sequence). Nested cycles are propagated in one run, other cycles step by
% step.

% System ------------------------
Sys.S = 1/2;
Sys.T1 = 0.2;
Sys.T2 = 0.1;

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.tp = 0.01;
Pulse.Flip = pi/2;

Exp.Sequence = {Pulse 0.1 Pulse};
Exp.mwFreq = 9.5;
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = [0 0.05];  % µs
Exp.DetTimeStep = 0.005;
Exp.DetOperator = {'z1','+1'};

% Options ---------------------------
Opt.SimFreq = 9.45;
Opt.IntTimeStep = 0.001;
Opt.Relaxation = 1;

tol = 1e-10;

% (a) cycle with zero mean weight on pulse 1
w = [1 1i -1 -1i];
Expa = Exp;
Expa.PhaseCycle = [(0:3).'*pi/2 zeros(4,1) -angle(w.')];
[~,signalPC] = saffron(Sys,Expa,Opt);

Expa = Exp;
Expa.nPoints = 4;
Expa.Dim1 = {'p1.Phase',pi/2};
[~,signalDim] = saffron(Sys,Expa,Opt);
signalRef = squeeze(sum(w.'.*signalDim,1))/4;

ok(1) = areequal(signalPC,signalRef,tol,'rel');

% (b) cycles with nonzero mean weights on pulses 1 and 2
w = [1 1i];
Expb = Exp;
[ph1,ph2] = ndgrid([0 pi/2]);
Expb.PhaseCycle = [ph1(:) ph2(:) -ph1(:)-ph2(:)];
[~,signalPC] = saffron(Sys,Expb,Opt);

Expb = Exp;
Expb.nPoints = [2 2];
Expb.Dim1 = {'p1.Phase',pi/2};
Expb.Dim2 = {'p2.Phase',pi/2};
[~,signalDim] = saffron(Sys,Expb,Opt);
signalRef = squeeze(sum(w.'.*w.*signalDim,[1 2]))/4;

ok(2) = areequal(signalPC,signalRef,tol,'rel');

% (c) both pulses stepped together (not a nested cycle), selecting the
% total coherence order change -1 (detected by S+)
w = [1 -1i -1 1i];
ph = (0:3).'*pi/2;
Expc = Exp;
Expc.PhaseCycle = [ph ph -angle(w.')];
[~,signalPC] = saffron(Sys,Expc,Opt);

Expc = Exp;
Expc.nPoints = 4;
Expc.Dim1 = {'p1.Phase,p2.Phase',pi/2};
[~,signalDim] = saffron(Sys,Expc,Opt);
signalRef = squeeze(sum(w.'.*signalDim,1))/4;

ok(3) = areequal(signalPC,signalRef,tol,'rel') && norm(signalRef(:,2))>0.1;

% (d) cycles (a) and (c) with detection over the whole sequence: the
% signal before and during the cycled pulses is phase-cycled, too (the
% zero mean weight of (a) removes the Sz signal before the first pulse)
Expd = Exp;
Expd.DetWindow = 'sequence';
ph = (0:3).'*pi/2;
cycles = {[1 1i -1 -1i] 'p1.Phase' [ph 0*ph]; [1 -1i -1 1i] 'p1.Phase,p2.Phase' [ph ph]};
for c = 1:2
  w = cycles{c,1};
  E = Expd;
  E.PhaseCycle = [cycles{c,3} -angle(w.')];
  [~,signalPC] = saffron(Sys,E,Opt);
  E = Expd;
  E.nPoints = 4;
  E.Dim1 = {cycles{c,2},pi/2};
  [~,signalDim] = saffron(Sys,E,Opt);
  signalRef = squeeze(sum(w.'.*signalDim,1))/4;
  ok(end+1) = areequal(signalPC,signalRef,tol,'abs') && max(abs(signalRef(:)))>0.1; %#ok<AGROW>
  if c==1
    ok(end+1) = abs(signalPC(1,1))<1e-12 && abs(signalDim(1,1,1))>0.1; %#ok<AGROW>
  end
end

warning(orig_state);
