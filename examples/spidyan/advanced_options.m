clear Exp Sys Opt Pulse
% Default Spin System
DefSys.S = 1/2;
DefSys.ZeemanFreq = 33.500; % GHz
DefSys.T1 = 1; % us
DefSys.T2 = 0.5; % us

% Pulse Definitions
Rectangular.Type = 'rectangular';

Adiabatic.Type = 'quartersin/linear';
Adiabatic.trise = 0.05; % us
Adiabatic.Qcrit = 5;

% A default Experiment/Sequence
DefExp.Field = 1240; % mT 
DefExp.TimeStep = 0.0001; % us
DefExp.mwFreq = 33.5; % GHz
DefExp.DetEvents = 1;

% Options
DefOpt.DetOperator = {'z1'};
DefOpt.FrameShift = 32; % GHz

%% 1) Using a custom initial state and equilibrium state:
% The spin relaxes to from the provided initial state to the equilibrium

Sys = DefSys;
Sys.initState = +sop(Sys.S,'z');
Sys.eqState = -sop(Sys.S,'z');

Exp = DefExp;
Exp.t = [0.2 5]; % us
Exp.Pulses = {Rectangular};
Exp.Frequency = 0;
Exp.Flip = pi/2;

Opt = DefOpt;
% Relaxation is active only during the second event
Opt.Relaxation = [0 1];

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% plotting
figure(1)
clf
plot(TimeAxis*1000,real(Signal))
xlabel('t [ns]')
ylabel('<S_i>')
legend(DefOpt.DetOperator)

%% 2) complex and custom excitation operators
% Complex excitation operators allows to investigate Bloch Siegert shifts
% If ComplexExcitation is active, the excitation operator takes the form:
% real(IQ)*Sx + imag(IQ)*Sy

Sys = DefSys;

Exp = DefExp;
Exp.t = 0.2; % us
Exp.Pulses = {Adiabatic};
Exp.Frequency = [-0.05 0.05]; % GHz

Opt = DefOpt;

[TimeAxis, Regular] = spidyan(Sys,Exp,Opt);

Opt.ComplexExcitation = 1;
[~, ComplexExOp] = spidyan(Sys,Exp,Opt);

% plotting
figure(2)
clf
hold on
plot(TimeAxis*1000,Regular)
plot(TimeAxis*1000,ComplexExOp)
xlabel('t [ns]')
ylabel('<S_i>')
legend('Regular ExOp','Complex ExOp')

% With custom excitation operators it is possible to investigate the effect
% of pulses during a pulse sequence in more detail - in this example our
% custom excitation operator is the same as the default one (Sx), but 
% manually defined (sop(Sys.S,'x') would also be possible. 
% Feel free to experiment!

Opt = DefOpt;
Opt.ExcOperator = {sop(Sys.S,'x(1|2)')};

[TimeAxis, Custom] = spidyan(Sys,Exp,Opt);

plot(TimeAxis*1000,Custom)
legend('Regular ExOp','Complex ExOp','Custom ExOp')

% It is also possible to combine ComplexExcitation with custom excitation 
% operators
% In this case the excitation operator takes the form:
% real(IQ)*real(ExOperator) + imag(IQ)*imag(ExOperator)

Opt = DefOpt;
Opt.ExcOperator = {sop(Sys.S,'x(1|2)')+sop(Sys.S,'y(1|2)')};
Opt.ComplexExcitation = 1;

[TimeAxis, Custom] = spidyan(Sys,Exp,Opt);

plot(TimeAxis*1000,Custom)
legend('Regular ExOp','Complex ExOp','Custom ExOp','Complex custom ExOp')

%% 3) state trajectories
% if state trajectories is switched on for an event, density matrices at
% each propagation point are stored and returned

Sys = DefSys;
Exp = DefExp;
Opt = DefOpt;

% Switches on state trajectories for both events
Opt.StateTrajectories = [1 1];

Exp.t = [0.2 0.5]; % us
Exp.Pulses = {Adiabatic};
Exp.Frequency = [-0.05 0.05]; % GHz

Exp.nPoints = [1 1];
Exp.Dim1 = {'p1.tp', 0.01};
Exp.Dim2 = {'p1.tp', 0.01};

[TimeAxis, Signal, AdvancedOutputs] = spidyan(Sys,Exp,Opt);

message = ['Cell array ''AdvancedOutputs.StateTrajectories'' contains ' num2str(size(AdvancedOutputs.StateTrajectories,2)) ' density matrices.'];
disp(message)
