clear Sys Exp Vary Opt Pulse sigmas Det

% System ------------------------
Sys.S = [1/2];
Sys.g = [1.9302];
Sys.ZeemanFreq = [33.500];  % New Field for Providing resonace frequeny(ies)
% Sys.initState = [0 0.5; 0.5 0]; % Field for the initial state, must recognize string or a matrix, collision with Exp.T?
% Sys.eqState = []; % Field for equilibrium state, same requirements as initial state
% Sys = nucspinadd(Sys,'14N',[14 14 32]);
% Sys.Nucs = '63Cu';
% Sys.A = [50 50 460 50 50 460];
% Sys.J = 0; 
Sys.T1 = 0.5;
Sys.T2 = 0.5;


% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
% Pulse.Qcrit = 0.5;
% Pulse.Amplitude = 30;
% Pulse.ComplexExcitation = true;

% PC = [0 1; pi -1];

Exp.t = [0.1 0.01 0.1 0.5 0.1 0.5];
Exp.Pulses = {Pulse 0 Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 1 1 1 1 1]; 
% Exp.PhaseCycle{1} = PC;

Exp.Resonator.nu0 = 33.5;
Exp.Resonator.QL = 200;
Exp.Resonator.Mode = 'compensate';
% Exp.Resonator.nu = [];
% Exp.Resonator.TransferFunction = [];

Exp.Dim1 = {'p2.Flip' -0.5};
Exp.Dim2 = {'p3.Position' -0.45};
% Exp.Dim2 = {'d2' 0.01};
Exp.nPoints = [2 2];

% Detection -------------------------
Opt.DetOperator = {'z1'};
Opt.FreqTranslation = [0 -33.5 -33.5]; 


% Options ---------------------------
Opt.Relaxation = [0 0];
% Opt.ExcOperator = {[0 0.5; 0.5 0] 'x1'};
Opt.FrameShift = 32;
Opt.Resonator = [];
Opt.StateTrajectories = [];


% Function Call -----------------------------

[t, signal, state, sigmas, Eventsnew] = spidyan(Sys,Exp,Opt);


% Plotting ----------------------------------
%% 
try
  figure(2);clf
  hold on
  for i = 1 : size(signal,1)
    plotsignal = squeeze(signal(i,:,:));
    plot(t(i,:),real(plotsignal))
  end
%   plot(t(i,:),imag(plotsignal))

catch
  figure(2);clf
  hold on
  for i = 1: length(signal)
    plot(t{i},real(signal{i}))
    
  end
end
ylabel('<{S_z}>')
xlabel('t [\mus]')