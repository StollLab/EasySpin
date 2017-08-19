clear Sys Exp Vary Opt Pulse sigmas 

Sys.S = 1/2;
% Sys.g = [2 2 2];
Sys.ResonanceFrequency = 1500;  % New Field for Providing resonace frequeny(ies)
Sys.InitialState = []; % Field for the initial state, must recognize string or a matrix, collision with Exp.T?
Sys.EquilibriumState = []; % Field for equilibrium state, same requirements as initial state
% Sys = nucspinadd(Sys,'14N',[14 14 32]);
Sys.T1 = 5;
Sys.T2 = 1;

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.ComplexExcitation = 0;
% PC = [0 1; pi -1];
PC = 0;

% Exp.t = [0.1 0.5 0.1 0.5 0.1 0.5 0.1];
% Exp.Pulses = {Pulse 0 Pulse 0 Pulse 0 Pulse};
Exp.t = [0.1 5];
Exp.Pulses = {Pulse};
Exp.B = 1240; % New Field: Magnetic Field

Exp.Frequency = [-100 100] + 1500;
Exp.Flip = [pi pi pi pi];
% Exp.PhaseCycle{1} = PC;
Exp.TimeStep = 0.0001; % us


% Opt.DetectionOperators = {'z1' 'x1' 'p1' 'm1' 'x2'}; % Need a field name here, make a new branch
% Opt.DetectionOperators = {'x1' 'p1'}; % Need a field name here, make a new branch
% Opt.DownConversionFrequency = [0 -1.5 -1.5 1.5]; 
Opt.DetectionOperators = {'z1'};
% Opt.ExcitationOperators = {[0 1; 0 0]}; %%%% Instead of using excitation
% operators, maybe move this field into pulse???

Opt.DetectedEvents = [1 1 1 1 1 1 1]; 
Opt.Relaxation = [0 1 0 0 0 0 0 0];
Opt.StateTrajectories = [];




% Exp.Inc1 = {'p2.Position,p3.Position' -0.2;
%             'd1' 0.2};
%           Exp.Inc2 = {'p2.Position,p3.Position' -0.3};
% Exp.nPoints = [4 2];

% Exp.Dim1 = {'p1.trise' .01};
% Exp.Dim1 = {'p3.Flip' -pi/4;
% Exp.Dim1 = {'p3.tp' 0.1 };
% Exp.Dim1 = {'d2' 0.2};
% 
% Exp.Dim1 = {'p3.Position' -0.35};
% Exp.Dim1 = {'p3.Position,p2.Position' -0.1};
%             'd1' .05};
% Exp.Dim2 = {'p2.Position' -0.2};
%             'd1',0.1};
% Exp.Dim2 = {'d1' 0.2};
% Exp.Dim3 = {'d1' 0};
% Exp.nPoints = [2];% 3];

[t, signal, state, sigmas, Eventsnew]=spidyan(Sys,Exp,Opt);
options.det_op = {[0 1; 0 0]};
options.awg.s_rate = 10;
options.dc = 1;

% [signal1] = strike(real(squeeze(signal))', t*1000, 1.5,options);
% [t,signal] = rfmixer(t,real(squeeze(signal)),1.5,'LSB');

% figure(1); clf
% plot(t,signal1)
% hold on
% plot(t,signal)

% disp(state)
% try
%   figure(1)
%   plot(t, squeeze(real(signal)))
% end

% whos signal
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