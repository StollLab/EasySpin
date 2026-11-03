function ok = test()

orig_state = warning;
warning('off','all');

% Phase cycling with relaxation: the cycled signal must equal the weighted
% sum of separate experiments with the pulse phases of the cycle steps,
% including the recovery toward equilibrium after the cycled pulses

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 9.5;
Sys.T1 = 0.2;
Sys.T2 = 0.1;

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.tp = 0.01;
Pulse.Flip = pi/2;

Exp.Sequence = {Pulse 0.1 Pulse 0.05};
Exp.mwFreq = 9.5;
Exp.DetSequence = [0 0 0 1];
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
Expa.PhaseCycle = {[(0:3).'*pi/2 w.']};
[~,signalPC] = spidyan(Sys,Expa,Opt);

Expa = Exp;
Expa.nPoints = 4;
Expa.Dim1 = {'p1.Phase',pi/2};
[~,signalDim] = spidyan(Sys,Expa,Opt);
signalRef = squeeze(sum(w.'.*signalDim,1))/4;

ok(1) = areequal(signalPC,signalRef,tol,'rel');

% (b) cycles with nonzero mean weights on pulses 1 and 2
w = [1 1i];
Expb = Exp;
Expb.PhaseCycle = {[0 1; pi/2 1i] [0 1; pi/2 1i]};
[~,signalPC] = spidyan(Sys,Expb,Opt);

Expb = Exp;
Expb.nPoints = [2 2];
Expb.Dim1 = {'p1.Phase',pi/2};
Expb.Dim2 = {'p2.Phase',pi/2};
[~,signalDim] = spidyan(Sys,Expb,Opt);
signalRef = squeeze(sum(w.'.*w.*signalDim,[1 2]))/4;

ok(2) = areequal(signalPC,signalRef,tol,'rel');

warning(orig_state);
