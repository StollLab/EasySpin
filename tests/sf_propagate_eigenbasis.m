function ok = test() %#ok<FNDEF>

% Free evolution is evaluated in the eigenbasis of the Hamiltonian, unless
% the density matrices are stored (Events{i}.StateTrajectories), which
% propagates step by step. Both agree, also for detection operators with
% several electron coherence orders: S = 1 with zero-field splitting, and
% two coupled electrons. FID after a pi/2 pulse, shifted simulation frame.

tol = 1e-9;

Exp.mwFreq = 9.5;  % GHz
Exp.Field = mhz2mt(9505);  % mT, resonance field for g = gfree at 9.505 GHz
p90.Flip = pi/2;
p90.tp = 0.02;  % µs
Exp.Sequence = {p90};
Exp.DetWindow = [0 0.1];  % µs
Exp.DetTimeStep = 1e-3;  % µs
Opt.IntTimeStep = 1e-5;  % µs
Opt.SimFreq = 9;  % GHz
Opt.DetFreq = Exp.mwFreq;

SysT.S = 1;
SysT.D = 30;  % MHz
OptT = Opt;
OptT.DetOperator = {'x1','+1','z1'};

SysE.S = [1/2 1/2];
SysE.g = gfree*[9.505; 9.51]/9.505;  % Zeeman frequencies 9.505 and 9.51 GHz
SysE.J = 5;  % MHz
OptE = Opt;
OptE.DetOperator = {'x1','x2','+1','+2','z1'};

[sT,sTstep] = signals(SysT,Exp,OptT);
[sE,sEstep] = signals(SysE,Exp,OptE);

ok(1) = areequal(sTstep,sT,tol,'abs');
ok(2) = areequal(sEstep,sE,tol,'abs');

end

% Signals from sf_propagate with eigenbasis and with step-by-step free
% evolution, set up as in the full-space engine of saffron (single
% orientation, molecular frame = lab frame)
function [sig,sigStep] = signals(Sys,Exp,Opt)
[Events,Vary,Opt] = runprivate('sf_buildevents',Exp,Opt);
[Sys,Sigma,Det,Events,Relaxation] = runprivate('sf_propagationsetup',Sys,Events,Opt);
H = ham(Sys,Exp.Field*[0 0 1]);
H = runprivate('sf_secularham',H,Sys.Spins,Sys.nElectrons);
Sz_total = 0;
for iSpin = 1:Sys.nElectrons
  Sz_total = Sz_total + sop(Sys.Spins,[iSpin,3]);
end
H = H - Opt.FrameShift*1e3*Sz_total;  % simulation frame (MHz)
% A single propagation run: no phase cycle (the only setting of each pulse)
% and no pathway selection
nPulses = sum(cellfun(@(ev)strcmp(ev.type,'pulse'),Events));
Runs.rows = repmat({1},1,nPulses);
Runs.factors = repmat({1},1,nPulses);
Runs.Projection = cell(1,nPulses);
Runs.weight = 1;
[~,sig] = runprivate('sf_propagate',Sigma,H,Det,Events,Relaxation,Vary,Runs);
for k = 1:numel(Events)
  Events{k}.StateTrajectories = true;
end
[~,sigStep] = runprivate('sf_propagate',Sigma,H,Det,Events,Relaxation,Vary,Runs);
end
