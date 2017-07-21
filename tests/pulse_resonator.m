function [err,data] = test(opt,olddata)
% Simulate the effect of the resonator on a rectangular pulse
%--------------------------------------------------------------------------

% Resonator defined through Q and resonance frequency
%--------------------------------------------------------------------------
Par.tp = 0.200; % us
Par.Type = 'rectangular';
Par.TimeStep = 0.00025;

Par.mwFreq = 9.5;

Opt.Resonator = 'simulation';
Par.ResonatorFrequency = Par.mwFreq;
Opt.CutoffFactor = 1/1000;

QL = 50:50:400;

for iq = 1:2:2*numel(QL)
  
  Par.QL = QL(round(iq/2));
  
  [t,IQ] = pulse(Par,Opt);
  
  % Fall time
  [~,ind] = min(abs(t-Par.tp));
  [k,c] = exponfit(t(ind:end),real(IQ(ind:end)),1,'noconst');
  tau_fall = 1/k; % tau = Q/(pi*nu_mw)
  Q_fall = tau_fall*pi*Par.ResonatorFrequency*1e3;
  
  suberr(iq) = (Par.QL-Q_fall)/Par.QL>0.025;
  
  % Rise time
  ind = numel(t)-ind;
  [k,c] = exponfit(t(1:ind),real(IQ(1:ind)),1);
  tau_rise = 1/k;
  Q_rise = tau_rise*pi*Par.ResonatorFrequency*1e3;
  
  suberr(iq+1) = (Par.QL-Q_rise)/Par.QL>0.025;
  
end

err(1) = any(suberr);

% Resonator profile given as input
%--------------------------------------------------------------------------
clearvars -except err

Par.tp = 0.200; % us
Par.Type = 'rectangular';
Par.TimeStep = 0.00025;

Par.mwFreq = 9.5;
Opt.Resonator = 'simulation';

QLvalues = 50:50:400;

for iq = 1:2:2*numel(QLvalues)
  
  QL = QLvalues(round(iq/2));
  Par.faxis = 9.2:0.010:9.8; % GHz
  Par.MagnitudeResponse = abs(1./(1+1i*QL*(Par.faxis/Par.mwFreq-Par.mwFreq./Par.faxis)));
  
  [t,IQ] = pulse(Par,Opt);
  
  % Fall time
  [~,ind] = min(abs(t-Par.tp));
  [k,c] = exponfit(t(ind:end),real(IQ(ind:end)),1,'noconst');
  tau_fall = 1/k; % tau = Q/(pi*nu_mw)
  Q_fall = tau_fall*pi*Par.mwFreq*1e3;
  
  suberr(iq) = (QL-Q_fall)/QL>0.05;
  
  % Rise time
  ind = numel(t)-ind;
  [k,c] = exponfit(t(1:ind),real(IQ(1:ind)),1);
  tau_rise = 1/k;
  Q_rise = tau_rise*pi*Par.mwFreq*1e3;
  
  suberr(iq+1) = (QL-Q_rise)/QL>0.05;
  
end

err(2) = any(suberr);

err = any(err);

data = [];

