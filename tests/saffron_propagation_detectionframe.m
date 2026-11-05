function ok = test() %#ok<FNDEF>

% Signals are returned in the detection frame (Exp.DetFreq, Exp.DetPhase),
% obtained exactly from the density matrix. Signals of detection operator D
% are Tr(D*rho)/Tr(D'*D). Single orientation, FID after a pi/2 pulse.

orig_state = warning;
warning('off','all');

Sys.S = 1/2;
Exp.mwFreq = 9.5;  % GHz
Exp.Field = mhz2mt(9505);  % mT, resonance field for g = gfree at 9.505 GHz
Exp.MolFrame = [0 0 0];  % single orientation
p90.Flip = pi/2;
p90.tp = 0.02;  % µs
Exp.Sequence = {p90};
Exp.DetWindow = [0 0.1];  % µs, FID after the pulse
Exp.DetTimeStep = 1e-3;  % µs
Exp.DetOperator = {'+1','-1','x1','y1','z1'};
Opt.IntTimeStep = 1e-5;  % µs
Opt.SimFreq = 9;

[t,s] = saffron(Sys,Exp,Opt);
plus = s(:,1); minus = s(:,2); x = s(:,3); y = s(:,4); z = s(:,5);

% (1) independent of the simulation frame, for both polarizations
Opt8 = Opt; Opt8.SimFreq = 8;
[~,s8] = saffron(Sys,Exp,Opt8);
OptOnRes = Opt; OptOnRes.SimFreq = Exp.mwFreq;
[~,sOnRes] = saffron(Sys,Exp,OptOnRes);
ExpC = Exp; ExpC.mwPolarization = 'circular';
[~,sC] = saffron(Sys,ExpC,Opt8);
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
[~,s2] = saffron(Sys,Exp2,Opt);
ok(3) = areequal(s2(:,5),z,1e-12,'abs');

% (4) DetPhase rotates the detection reference
phi = pi/3;
Exp3 = Exp;
Exp3.DetPhase = phi;
[~,s3] = saffron(Sys,Exp3,Opt);
ok(4) = areequal(s3(:,1),exp(-1i*phi)*plus,1e-12,'abs') && ...
  areequal(s3(:,2),exp(1i*phi)*minus,1e-12,'abs') && ...
  areequal(s3(:,3),cos(phi)*x+sin(phi)*y,1e-12,'abs');

% (5) the same phase added to all pulses and to the detection does not
% change the signal (initial state without electron coherences)
alpha = 0.7;
Exp4 = Exp;
Exp4.Sequence{1}.Phase = alpha;
Exp4.DetPhase = alpha;
[~,s4] = saffron(Sys,Exp4,Opt);
ok(5) = areequal(s4,s,1e-10,'abs');

% (6) DetFreq = mwFreq + delta multiplies S+ by exp(-2i*pi*delta*t), with
% t from the start of the sequence (the time axis of a detection window
% starts at the end of the sequence)
delta = 0.01;  % GHz
Exp5 = Exp;
Exp5.DetFreq = Exp.mwFreq + delta;
[~,s5] = saffron(Sys,Exp5,Opt);
ok(6) = areequal(s5(:,1),plus.*exp(-2i*pi*delta*1e3*(t(:)+p90.tp)),1e-12,'abs');

% (7) a single point equals the transient at the same time
Exp6 = Exp;
Exp6.Sequence = {p90 0.05};
Exp6.DetWindow = 0;
s6 = saffron(Sys,Exp6,Opt);
[~,idx] = min(abs(t-0.05));
ok(7) = areequal(s6(:).',s(idx,:),1e-9,'abs');

% (8) the single point at the end of the sequence does not depend on the
% simulation frame, and equals the last point of the transient
Exp7 = Exp;
Exp7.Sequence = {p90 0.1};
Exp7.DetWindow = 0;
sEnd = saffron(Sys,Exp7,Opt);
sEnd8 = saffron(Sys,Exp7,Opt8);
ok(8) = areequal(sEnd8,sEnd,1e-3,'abs') && areequal(sEnd(:).',s(end,:),1e-9,'abs');

% (9) invalid frames
ExpNeg = Exp; ExpNeg.DetFreq = -9.5;
ExpVec = Exp; ExpVec.DetOperator = {'z1','+1'}; ExpVec.DetFreq = [0 9.5];
OptNeg = Opt; OptNeg.SimFreq = -1;
ok(9) = throwsWith(@()saffron(Sys,ExpNeg,Opt),'DetFreq') && ...
  throwsWith(@()saffron(Sys,ExpVec,Opt),'single frequency') && ...
  throwsWith(@()saffron(Sys,Exp,OptNeg),'SimFreq');

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
