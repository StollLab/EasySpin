function ok = test()

% Sys.initState in all input forms, Hilbert-space propagation

tol = 1e-10;

Pulse.Type = 'rectangular';
Pulse.tp = 0.01;
Pulse.Flip = pi/2;

Exp.Sequence = {Pulse 0.05};
Exp.Field = 340;
Exp.mwFreq = 9.5;
Exp.DetSequence = [0 1];
Opt.IntTimeStep = 5e-4;
Opt.SimFreq = 9.48;

% Fully polarized S=1/2, compared to the default initial state -Sz
Sys.S = 1/2;
Sys.ZeemanFreq = 9.5;
[~,ref] = spidyan(Sys,Exp,Opt);
Sys.initState = {[1 0],'eigen'};
[~,sig1] = spidyan(Sys,Exp,Opt);
Sys.initState = {[0 1],'uncoupled'};
[~,sig2] = spidyan(Sys,Exp,Opt);

% Electron-only state with a nucleus: trace is kept, so -Sz/nNuc
SysN.S = 1/2;
SysN.g = 2;
SysN.Nucs = '1H';
SysN.A = [5 10];
ExpN = Exp;
ExpN.Field = 339.4;
[~,refN] = spidyan(SysN,ExpN,Opt);
SysN.initState = {diag([0 1]),'uncoupled'};
[~,sigN] = spidyan(SysN,ExpN,Opt);

% Triplet: 'xyz' populations equal 'zerofield' populations in energy order
% (D>0, E>0: Z < X < Y)
SysT.S = 1;
SysT.ZeemanFreq = 9.5;
SysT.D = [300 30];
SysT.initState = {[0.5 0.3 0.2],'xyz'};
[~,sigT1] = spidyan(SysT,Exp,Opt);
SysT.initState = {[0.2 0.5 0.3],'zerofield'};
[~,sigT2] = spidyan(SysT,Exp,Opt);

% Populations and coherences at t=0 via detection operators
ExpD = Exp;
ExpD.DetSequence = [1 0];
ExpD.DetOperator = {'a1','b1'};
Sys.initState = {[0.2 0.8],'eigen'};
[~,sigP] = spidyan(Sys,ExpD,Opt);
ExpD.DetOperator = {'x1','y1'};
Sys.initState = sop(1/2,'x');
[~,sigC] = spidyan(Sys,ExpD,Opt);

ok(1) = areequal(sig1,ref,tol,'rel');
ok(2) = areequal(sig2,ref,tol,'rel');
ok(3) = areequal(sigN,refN/2,tol,'rel');
ok(4) = areequal(sigT1,sigT2,tol,'rel');
ok(5) = areequal(sigP(1,:),[0.8 0.2],tol,'abs');
ok(6) = areequal(sigC(1,:),[1 0],tol,'abs');

