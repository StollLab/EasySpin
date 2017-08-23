clear Sys Exp Vary Opt Pulse sigmas 

% System ------------------------
Sys.S = [1];
Sys.g = [2 2 2];
Sys.ZeemanFreq = [1500];  % New Field for Providing resonace frequeny(ies)
% Sys.initState = [0 0.5; 0.5 0]; % Field for the initial state, must recognize string or a matrix, collision with Exp.T?
% Sys.eqState = []; % Field for equilibrium state, same requirements as initial state
% Sys = nucspinadd(Sys,'14N',[14 14 32]);
% Sys.Nucs = '63Cu';
% Sys.A = [50 50 460 50 50 460];
% Sys.J = 0; 
Sys.T1 = 0.5;
% Sys.T2 = 0.5;


% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.ComplexExcitation = 0;

% PC = [0 1; pi -1];

Exp.t = [0.1 0.5 0.1 0.5 0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse 0 Pulse 0 Pulse};
Exp.B = 1240; % New Field: Magnetic Field
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-100 100] + 1500;
Exp.Flip = [pi pi pi pi];
% Exp.PhaseCycle{1} = PC;

% Exp.Dim1 = {'p2.Position' -0.2};
%             'd1',0.1};
% Exp.nPoints = [2];

% Options ---------------------------
Opt.DetectionOperators = {'z1' 'x1' 'p1' 'm1'}; % Need a field name here, make a new branch
Opt.FreqTranslation = [0 -1.5 -1.5 1.5]; 
Opt.Relaxation = [0 0 0 0 0 0 0 0];
% Opt.ExcitationOperators = {[0 1; 0 0]};

% Move this into detection structure?
Opt.DetectedEvents = [1 1 1 1 1 1 1]; 
Opt.StateTrajectories = [];


% Function Call -----------------------------

[t, signal, state, sigmas, Eventsnew]=spidyan(Sys,Exp,Opt);


% Plotting ----------------------------------

try
  figure(2);clf
  hold on
  for i = 1 : size(signal,1)
    plotsignal = squeeze(signal(i,:,:));
    plot(t(i,:),(plotsignal))
  end

catch
  figure(2);clf
  hold on
  for i = 1: length(signal)
    plot(t{i},real(signal{i}))
    
  end
end
ylabel('<{S_z}>')
xlabel('t [\mus]')