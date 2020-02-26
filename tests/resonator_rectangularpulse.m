function ok = test()

% Check the simulation of the effect of the resonator on a rectangular
% pulse by comparing the Q values estimated from the rise and fall
% times of the simulated pulse with the Q values used for the resonator
% definition
%--------------------------------------------------------------------------

% Resonator defined through Q and resonance frequency
%--------------------------------------------------------------------------
Par.tp = 0.200; % us
Par.Type = 'rectangular';
Par.TimeStep = 0.00025;

[t,IQ] = pulse(Par);

mwFreq = 9.5;
QL = 50:50:400;

for iq = 1:2:2*numel(QL)
  
  [t_,IQ_] = resonator(t,IQ,mwFreq,mwFreq,QL(round(iq/2)),'simulate');
  
  % Fall time
  [ignore,ind] = min(abs(t_-Par.tp));
  [k,c] = exponfit(t_(ind:end),real(IQ_(ind:end)),1,'noconst');
  tau_fall = 1/k; % tau = Q/(pi*nu_mw)
  Q_fall = tau_fall*pi*mwFreq*1e3;
  
  suberr(iq) = (QL(round(iq/2))-Q_fall)/QL(round(iq/2))>0.025;
  
  % Rise time
  ind = numel(t_)-ind;
  [k,c] = exponfit(t_(1:ind),real(IQ_(1:ind)),1);
  tau_rise = 1/k;
  Q_rise = tau_rise*pi*mwFreq*1e3;
  
  suberr(iq+1) = (QL(round(iq/2))-Q_rise)/QL(round(iq/2))>0.025;
  
end

ok(1) = ~any(suberr);

% Resonator profile given as input
%--------------------------------------------------------------------------
clear Par

Par.tp = 0.200; % us
Par.Type = 'rectangular';
Par.TimeStep = 0.00025;

[t,IQ] = pulse(Par);

mwFreq = 9.5;

QLvalues = 50:50:400;

for iq = 1:2:2*numel(QLvalues)
  
  QL = QLvalues(round(iq/2));
  f = 9.2:0.010:9.8; % GHz
  H = abs(1./(1+1i*QL*(f/mwFreq-mwFreq./f)));
  
  [t_,IQ_] = resonator(t,IQ,mwFreq,f,H,'simulate');
  
  % Fall time
  [ignore,ind] = min(abs(t_-Par.tp));
  [k,c] = exponfit(t_(ind:end),real(IQ_(ind:end)),1,'noconst');
  tau_fall = 1/k; % tau = Q/(pi*nu_mw)
  Q_fall = tau_fall*pi*mwFreq*1e3;
  
  suberr(iq) = (QL-Q_fall)/QL>0.05;
  
  % Rise time
  ind = numel(t_)-ind;
  [k,c] = exponfit(t_(1:ind),real(IQ_(1:ind)),1);
  tau_rise = 1/k;
  Q_rise = tau_rise*pi*mwFreq*1e3;
  
  suberr(iq+1) = (QL-Q_rise)/QL>0.05;
  
end

ok(2) = ~any(suberr);

% Transfer function given as input
%--------------------------------------------------------------------------
clear Par

Par.tp = 0.200; % us
Par.Type = 'rectangular';
Par.TimeStep = 0.00025;

[t,IQ] = pulse(Par);

mwFreq = 9.5;

QLvalues = 50:50:400;

for iq = 1:2:2*numel(QLvalues)
  
  QL = QLvalues(round(iq/2));
  f = 9.2:0.010:9.8; % GHz
  H = 1./(1+1i*QL*(f/mwFreq-mwFreq./f));
  
  [t_,IQ_] = resonator(t,IQ,mwFreq,f,H,'simulate');
  
  % Fall time
  [ignore,ind] = min(abs(t_-Par.tp));
  [k,c] = exponfit(t_(ind:end),real(IQ_(ind:end)),1,'noconst');
  tau_fall = 1/k; % tau = Q/(pi*nu_mw)
  Q_fall = tau_fall*pi*mwFreq*1e3;
  
  suberr(iq) = (QL-Q_fall)/QL>0.05;
  
  % Rise time
  ind = numel(t_)-ind;
  [k,c] = exponfit(t_(1:ind),real(IQ_(1:ind)),1);
  tau_rise = 1/k;
  Q_rise = tau_rise*pi*mwFreq*1e3;
  
  suberr(iq+1) = (QL-Q_rise)/QL>0.05;
  
end

ok(3) = ~any(suberr);
