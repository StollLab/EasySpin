function ok = test() %#ok<FNDEF>

% Sys.initState in Hilbert-space propagation, single orientation: a
% polarized state in the uncoupled basis, an electron-only state with a
% nucleus, and populations and coherences at t = 0. ('eigen', 'xyz' and
% 'zerofield' are tested in saffron_propagation_initstate.)

tol = 1e-10;

Pulse.Type = 'rectangular';
Pulse.tp = 0.01;
Pulse.Flip = pi/2;

Exp.Sequence = {Pulse};
Exp.Field = 340;
Exp.mwFreq = 9.5;
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = [0 0.05];  % µs
Opt.IntTimeStep = 5e-4;
Opt.SimFreq = 9.48;

% Fully polarized S=1/2, compared to the default initial state -Sz
Sys.S = 1/2;
Sys.g = gfree*mhz2mt(9500)/Exp.Field;  % Zeeman frequency 9.5 GHz
[~,ref] = saffron(Sys,Exp,Opt);
Sys.initState = {[0 1],'uncoupled'};
[~,sig] = saffron(Sys,Exp,Opt);

% Electron-only state with a nucleus: trace is kept, so -Sz/nNuc
SysN.S = 1/2;
SysN.g = 2;
SysN.Nucs = '1H';
SysN.A = [5 10];
ExpN = Exp;
ExpN.Field = 339.4;
[~,refN] = saffron(SysN,ExpN,Opt);
SysN.initState = {diag([0 1]),'uncoupled'};
[~,sigN] = saffron(SysN,ExpN,Opt);

% Populations and coherences at t=0 via detection operators
ExpD = Exp;
ExpD.DetWindow = 'sequence';
ExpD.DetOperator = {'a1','b1'};
Sys.initState = {[0.2 0.8],'eigen'};
[~,sigP] = saffron(Sys,ExpD,Opt);
ExpD.DetOperator = {'x1','y1'};
Sys.initState = sop(1/2,'x');
[~,sigC] = saffron(Sys,ExpD,Opt);

ok(1) = areequal(sig,ref,tol,'rel');
ok(2) = areequal(sigN,refN/2,tol,'rel');
ok(3) = areequal(sigP(1,:),[0.8 0.2],tol,'abs');
ok(4) = areequal(sigC(1,:),[1 0],tol,'abs');
