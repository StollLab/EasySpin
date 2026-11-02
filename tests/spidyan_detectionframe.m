function ok = test()

% Signals and density matrices are returned in the detection frame
% (Exp.DetFreq, Exp.DetPhase), obtained exactly from the density matrix.
% Signals of detection operator D are Tr(D*rho)/Tr(D'*D).

orig_state = warning;
warning('off','all');

Sys.S = 1/2;
Sys.ZeemanFreq = 9.505;  % GHz
Exp.mwFreq = 9.5;  % GHz
p90.Flip = pi/2;
p90.tp = 0.02;  % µs
Exp.Sequence = {p90 0.1};
Exp.DetSequence = [0 1];
Exp.DetTimeStep = 1e-3;  % µs
Exp.DetOperator = {'+1','-1','x1','y1','z1'};
Opt.IntTimeStep = 1e-5;  % µs
Opt.SimFreq = 9;

[t,s] = spidyan(Sys,Exp,Opt);
plus = s(:,1); minus = s(:,2); x = s(:,3); y = s(:,4); z = s(:,5);

% (1) independent of the simulation frame, for both polarizations
Opt8 = Opt; Opt8.SimFreq = 8;
[~,s8] = spidyan(Sys,Exp,Opt8);
OptOnRes = Opt; OptOnRes.SimFreq = Exp.mwFreq;
[~,sOnRes] = spidyan(Sys,Exp,OptOnRes);
ExpC = Exp; ExpC.mwPolarization = 'circular';
[~,sC] = spidyan(Sys,ExpC,Opt8);
ok(1) = areequal(s8,s,1e-3,'abs') && areequal(sOnRes,s,1e-3,'abs') && ...
  areequal(sC,s,1e-3,'abs');

% (2) <S-> = conj(<S+>); Sx and Sy are real (normalization Tr(D'*D):
% Sx = 2*Re(S+), Sy = 2*Im(S+) for S = 1/2)
ok(2) = areequal(minus,conj(plus),1e-12,'abs') && isreal(x) && isreal(y) && ...
  areequal(x,2*real(plus),1e-12,'abs') && areequal(y,2*imag(plus),1e-12,'abs');

% (3) Sz does not depend on DetFreq and DetPhase
Exp2 = Exp;
Exp2.DetFreq = 9.6;
Exp2.DetPhase = [1 2 3 4 5];
[~,s2] = spidyan(Sys,Exp2,Opt);
ok(3) = areequal(s2(:,5),z,1e-12,'abs');

% (4) DetPhase rotates the detection reference
phi = pi/3;
Exp3 = Exp;
Exp3.DetPhase = phi*ones(1,5);
[~,s3] = spidyan(Sys,Exp3,Opt);
ok(4) = areequal(s3(:,1),exp(-1i*phi)*plus,1e-12,'abs') && ...
  areequal(s3(:,2),exp(1i*phi)*minus,1e-12,'abs') && ...
  areequal(s3(:,3),cos(phi)*x+sin(phi)*y,1e-12,'abs');

% (5) the same phase added to all pulses and to the detection does not
% change the signal (initial state without electron coherences)
alpha = 0.7;
Exp4 = Exp;
Exp4.Sequence{1}.Phase = alpha;
Exp4.DetPhase = alpha*ones(1,5);
[~,s4] = spidyan(Sys,Exp4,Opt);
ok(5) = areequal(s4,s,1e-10,'abs');

% (6) DetFreq = mwFreq + delta multiplies S+ by exp(-2i*pi*delta*t)
delta = 0.01;  % GHz
Exp5 = Exp;
Exp5.DetFreq = Exp.mwFreq + delta;
[~,s5] = spidyan(Sys,Exp5,Opt);
ok(6) = areequal(s5(:,1),plus.*exp(-2i*pi*delta*1e3*t(:)),1e-12,'abs');

% (7) a single point equals the transient at the same time
Exp6 = Exp;
Exp6.Sequence = {p90 0.05 0};
Exp6.DetSequence = [0 0 1];
s6 = spidyan(Sys,Exp6,Opt);
[~,idx] = min(abs(t-(p90.tp+0.05)));
ok(7) = areequal(s6(:).',s(idx,:),1e-9,'abs');

% (8) the eigenbasis evaluation of free evolution and the step-by-step
% propagation (with state trajectories) agree, also for detection
% operators with several electron coherence orders
SysT.S = 1;
SysT.ZeemanFreq = 9.505;
SysT.D = 30;  % MHz
ExpT = Exp;
ExpT.DetOperator = {'x1','+1','z1'};
SysE.S = [1/2 1/2];
SysE.ZeemanFreq = [9.505 9.51];
SysE.J = 5;  % MHz
ExpE = Exp;
ExpE.DetOperator = {'x1','x2','+1','+2','z1'};
OptStep = Opt;
OptStep.StateTrajectories = true;
[~,sT] = spidyan(SysT,ExpT,Opt);
[~,sTstep] = spidyan(SysT,ExpT,OptStep);
[~,sE] = spidyan(SysE,ExpE,Opt);
[~,sEstep] = spidyan(SysE,ExpE,OptStep);
ok(8) = areequal(sT,sTstep,1e-9,'abs') && areequal(sE,sEstep,1e-9,'abs');

% (9) the final state is in the detection frame: independent of the
% simulation frame, and consistent with the last signal point
[~,~,info] = spidyan(Sys,Exp,Opt);
[~,~,info8] = spidyan(Sys,Exp,Opt8);
D = sop(Sys.S,'+1');
lastPoint = trace(D*info.FinalState)/trace(D'*D);
ok(9) = areequal(info8.FinalState,info.FinalState,1e-3,'abs') && ...
  abs(lastPoint-plus(end)) < 1e-12;

% (10) invalid frames
ExpNeg = Exp; ExpNeg.DetFreq = -9.5;
ExpVec = Exp; ExpVec.DetOperator = {'z1','+1'}; ExpVec.DetFreq = [0 9.5];
OptNeg = Opt; OptNeg.SimFreq = -1;
ok(10) = throwsWith(@()spidyan(Sys,ExpNeg,Opt),'DetFreq') && ...
  throwsWith(@()spidyan(Sys,ExpVec,Opt),'single frequency') && ...
  throwsWith(@()spidyan(Sys,Exp,OptNeg),'SimFreq');

warning(orig_state);

end

function rejected = throwsWith(fun,str)
try
  fun();
  rejected = false;
catch err
  rejected = contains(err.message,str);
end
end
