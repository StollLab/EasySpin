clear Sys Exp Vary Opt Pulse sigmas 
Sys = [];
% Exp = [];
Vary = [];
Opt = [];

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.ComplexExcitation = 0;

% PC = [0 1; pi -1];
PC = 0;
Exp.t = [0.1 0.5 0.1 0.5 0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse 0 Pulse 0 Pulse};


Opt.Detection = [1 1 1 1 1 1 1];
Opt.Relaxation = [0 0 0 0 0 0 0 0];
Opt.StateTrajectories = [];
Opt.DetectionOperator = [];
Opt.ExcitationOperator = [];

Exp.Frequency = [-100 100] +1500;
Exp.Flip = [pi pi pi pi];
Exp.PhaseCycle{1} = PC;

Exp.TimeStep = 0.0001; % us

% Exp.Inc1 = {'p2.Position,p3.Position' -0.2;
%             'd1' 0.2};
%           Exp.Inc2 = {'p2.Position,p3.Position' -0.3};
% Exp.nPoints = [4 2];

% Exp.Dim1 = {'d2' .2};
Exp.Dim1 = {'p3.Position' -0.35};
% Exp.Dim1 = {'p3.Position,p2.Position' -0.1};
%             'd1' .05};
% Exp.Dim2 = {'p2.Position' -0.2};
%             'd1',0.1};
% Exp.Dim2 = {'d1' 0.2};
% Exp.Dim3 = {'d1' 0};
Exp.nPoints = [4];% 3];

[t, signal, state, sigmas, Eventsnew]=spidyan(Sys,Exp,Opt);


% disp(state)
% try
%   figure(1)
%   plot(t, squeeze(real(signal)))
% end

whos signal
try
a = squeeze(signal(1,:,:));
b = squeeze(signal(2,:,:));
c = squeeze(signal(3,:,:));
% d = squeeze(signal(4,:,:));
figure(2); clf
plot(t(1,:),real(a))
hold on
plot(t(2,:),real(b))
plot(t(3,:),real(c))
% plot(t(4,:),real(d))

catch
  figure(4);clf
  hold on
  for i = 1: length(signal)
    plot(t{i},real(signal{i}))
    
  end
end

% 
% Events{1}.Detection = 0;
% Events{3}.Detection = 0;
% Events{2}.Detection = 0;
% Events{4}.Detection = 0;
%         
% tic
% [t, signal,state,sigma,Eventnew]=evolve2(Sigma,Ham, Det, Events, Relaxation);
% toc
% % profile on
% tic
% [t, signal,state,sigmas]=evolve2(Sigma,Ham, Det, Eventsnew, Relaxation, Vary);
% % profile off
% toc
% % profile report
% disp(state)
% % figure(2); clf
% % plot(t,real(signal))
% % 
% Events{1}.Detection = 0;
% Events{3}.Detection = 0;
% Events{2}.Detection = 0;
% Events{4}.Detection = 0;
%         
% tic
% [t, signal,state,sigma,Eventsnew]=evolve2(Sigma,Ham, Det, Events, Relaxation,Vary);
% toc
% disp(state)
% 
% tic
% [t, signal,state,sigma,Eventsnew]=evolve2(Sigma,Ham, Det, Eventsnew, Relaxation,Vary);
% toc
% disp(state)

ylabel('<{S_z}>')
xlabel('t [\mus]')