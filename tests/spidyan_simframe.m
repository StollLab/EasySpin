function ok = test()

% In a shifted simulation frame, pulses contain only the field components
% that rotate with the frame (rotating-wave approximation). The results
% must not depend on the frame frequency, also for linear polarization and
% for a frame rotating at the microwave frequency, and must agree with the
% lab frame up to counter-rotating effects (nu1/nu0). Excitation operators
% that a microwave pulse cannot drive in a shifted frame are rejected.

orig_state = warning;
warning('off','all');

Sys.S = 1/2;
Sys.ZeemanFreq = 9.505;  % GHz
Exp.Field = 340;
Exp.mwFreq = 9.5;

% Two pi/2 pulses with a varying relative phase, separated by a delay
Pulse.tp = 0.02;
Pulse.Flip = pi/2;
Exp.DetOperator = {'z1'};
phases = (0:3)*pi/2;
Sz = sop(Sys.S,'z');

% final electron polarization for each phase of the second pulse
function mz = finalSz(Exp,Opt)
  mz = zeros(size(phases));
  for k = 1:numel(phases)
    Pulse2 = Pulse;
    Pulse2.Phase = phases(k);
    Exp.Sequence = {Pulse 0.0123 Pulse2};
    [~,~,info] = spidyan(Sys,Exp,Opt);
    mz(k) = real(trace(Sz*info.FinalState));
  end
end

Opt.IntTimeStep = 1e-5;  % µs
Opt.SimFreq = 9;
mzLinear = finalSz(Exp,Opt);
Opt.SimFreq = 8;
mzLinear8 = finalSz(Exp,Opt);
Opt.SimFreq = Exp.mwFreq;
mzLinearOnRes = finalSz(Exp,Opt);
Opt.SimFreq = 9;
ExpC = Exp;
ExpC.mwPolarization = 'circular';
mzCircular = finalSz(ExpC,Opt);

Opt.SimFreq = 0;
Opt.IntTimeStep = 1e-6;  % µs
mzLab = finalSz(Exp,Opt);

ok(1) = areequal(mzLinear,mzLinear8,1e-3,'abs');
ok(2) = areequal(mzLinear,mzLinearOnRes,1e-3,'abs');
ok(3) = areequal(mzLinear,mzCircular,1e-12,'abs');
ok(4) = areequal(mzLinear,mzLab,2e-3,'abs');
% the phase of the second pulse matters
ok(5) = max(mzLinear)-min(mzLinear) > 0.5;

% Excitation operators with elements that do not change the electron mS by
% +-1 are rejected in a shifted frame
Exp.Sequence = {Pulse};
Opt.SimFreq = 9;
Opt.IntTimeStep = 1e-5;
SysN = Sys;
SysN.Nucs = '1H';
SysN.A = 5;  % MHz
OptN = Opt;
OptN.ExcOperator = {'x2'};
SysT = Sys;
SysT.S = 1;
OptT = Opt;
OptT.ExcOperator = {'x(1|3)'};
ok(6) = throwsError(@() spidyan(SysN,Exp,OptN));
ok(7) = throwsError(@() spidyan(SysT,Exp,OptT));

% An electron that precesses in the negative sense (negative g) is
% rejected in a shifted frame
SysG.S = [1/2 1/2];
SysG.ZeemanFreq = [9.505 -9.6];  % GHz
SysG.J = 1;  % MHz
ok(8) = throwsError(@() spidyan(SysG,Exp,Opt));

warning(orig_state);

end

function rejected = throwsError(fun)
try
  fun();
  rejected = false;
catch e
  rejected = contains(e.message,'shifted simulation frame');
end
end
