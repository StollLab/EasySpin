% Illustration of amplitude compression and compensation for transmitter
% nonlinearity
%==========================================================================
% In this example, two shaped pulses with different amplitudes are computed
% and the effect of transmitter nonlinearity is analyzed. The ideal pulse 
% inversion profiles are compared to the excitation profiles obtained for
% the predicted pulse output in the presence on amplitude compression.
% The pulse shapes are then compensated for amplitude distortions and the
% effect of the transmitter on the compensated pulses is modelled. The 
% output obtained for the compensated pulse shapes is in good agreement
% with the desired ideal pulses.

clear, clc, clf

% Ideal pulses
% ---------------------------------------------------------------------------------------
PulseAmplitudes = [13 60]; % MHz
delta = 0.300; % µs

% Pulse 1
Params1.Type = 'sech/tanh';
Params1.TimeStep = 0.00025; % µs
Params1.Frequency = [-50 50]; % MHz
Params1.tp = 0.200; % µs
Params1.beta = 10;
Params1.Amplitude = PulseAmplitudes(1); % MHz

[t1,signal1_ideal] = pulse(Params1);
[offset1_ideal,M1_ideal] = exciteprofile(t1,signal1_ideal);

% Pulse 2
Params2 = Params1;
Params2.tp = 0.100; % µs
Params2.Amplitude = PulseAmplitudes(2); % MHz
[t2,signal2_ideal] = pulse(Params2);
[offset2_ideal,M2_ideal] = exciteprofile(t2,signal2_ideal);

% Definition of transmitter nonlinearity
% ---------------------------------------------------------------------------------------
Ain =  0:0.001:1;
Gain = 85.5; % extrapolated output amplitude for max input amplitude in linear regime
Input1dBPoint = 0.60;
Aout = (Gain*Ain - (Gain*(1-10^(-1/20))/Input1dBPoint^2)*Ain.^3);
[~,ind] = min(abs(diff(Aout)));
Asat = Aout(ind);
Aout(ind:end) = Asat;

% Amplitude compression
% ---------------------------------------------------------------------------------------
% Pulse 1
signal1_compressed = transmitter(signal1_ideal/max(PulseAmplitudes),Ain,Aout,'simulate');

[offset1_compressed,M1_compressed] = exciteprofile(t1,signal1_compressed);

% Pulse 2
signal2_compressed = transmitter(signal2_ideal/max(PulseAmplitudes),Ain,Aout,'simulate');

[offset2_compressed,M2_compressed] = exciteprofile(t2,signal2_compressed);

% Amplitude compensation
% ---------------------------------------------------------------------------------------
nu1_out = max(Aout);

% Pulse 1
signal1_compensated = transmitter(signal1_ideal*nu1_out/max(PulseAmplitudes),Ain,Aout,'compensate');

% Pulse 2
signal2_compensated = transmitter(signal2_ideal*nu1_out/max(PulseAmplitudes),Ain,Aout,'compensate');

% Effect of the transmitter on compensated pulses
% ---------------------------------------------------------------------------------------

% Pulse 1
signal1_corrected = transmitter(signal1_compensated,Ain,Aout,'simulate');

[offset1,M1_corrected] = exciteprofile(t1,signal1_corrected);

% Pulse 2
signal2_corrected = transmitter(signal2_compensated,Ain,Aout,'simulate');

[offset2_corrected,M2_corrected] = exciteprofile(t2,signal2_corrected);

% Plot results
% ---------------------------------------------------------------------------------------
% Transmitter characterization
subplot(2,2,1)
hold on; box on; grid on;
title('Transmitter nonlinearity characterization')
plot(Ain,Ain*Gain,'--k')
plot(Ain,Aout,'b')
line([1 1]*Input1dBPoint,ylim,'LineStyle','--','Color',[0.7 0.7 0.7])
line(xlim,[1 1]*Gain*Input1dBPoint*10^(-1/20),'LineStyle','--','Color',[0.7 0.7 0.7])
text(0.02,Gain*Input1dBPoint*10^(-1/20)+5,'1 dB compression','FontSize',7)
line(xlim,[1 1]*max(Aout),'LineStyle','--','Color','k')
text(0.02,max(Aout)+5,'saturation','FontSize',7)
xlabel('Input amplitude (norm.)')
ylabel('Output amplitude (MHz)')

% Pulse shapes
subplot(2,2,2)
hold on; box on; grid on;
title('Pulses')
plot(t1,real(signal1_ideal),'Color',[0.7 0.7 0.7])
plot(t1,real(signal1_compressed),'r')
plot(t1,real(signal1_corrected),'b')
legend('ideal','compressed','corrected','Location','Best')
legend boxoff
plot(t2+delta,real(signal2_ideal),'Color',[0.7 0.7 0.7])
plot(t2+delta,real(signal2_compressed),'r')
plot(t2+delta,real(signal2_corrected),'b')
axis([0 0.4 -65 65])
line(xlim,[0 0],'Color','k')
line(xlim,[1 1]*max(Aout),'Color','k')
xlabel(['t (',char(181),'s)'])
ylabel('\nu_1 (MHz)')

% Compensated pulse shapes
subplot(2,2,3)
hold on; box on; grid on;
title('Compensated pulse shapes')
plot(t1,real(signal1_ideal)/max(real(signal2_ideal)),'Color',[0.7 0.7 0.7])
plot(t1,real(signal1_compensated),'b')
legend('ideal','compensated','Location','Best')
legend boxoff
plot(t2+delta,real(signal2_ideal)/max(real(signal2_ideal)),'Color',[0.7 0.7 0.7])
plot(t2+delta,real(signal2_compensated),'b')
axis([0 0.4 -1 1])
line(xlim,[0 0],'Color','k')
xlabel(['t (',char(181),'s)'])
ylabel('Input amplitude (norm.)')

% Excitation profiles
subplot(2,2,4)
hold on; box on; grid on;
title('Excitation profiles')
plot(offset1_ideal,M1_ideal(3,:),'Color',[0.7 0.7 0.7])
plot(offset1_compressed,M1_compressed(3,:),'r')
plot(offset1,M1_corrected(3,:),'b')
legend('ideal','compressed','corrected','Location','Best')
legend boxoff
plot(offset2_ideal,M2_ideal(3,:),'--','Color',[0.7 0.7 0.7])
plot(offset2_compressed,M2_compressed(3,:),'--r')
plot(offset2_corrected,M2_corrected(3,:),'--b')
xlabel('\Delta\nu (MHz)')
ylabel('M_z/M_0')
