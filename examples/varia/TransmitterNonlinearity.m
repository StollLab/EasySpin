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
delta = 0.300; % us

% Pulse 1
Params1.Type = 'sech/tanh';
Params1.TimeStep = 0.00025; % us
Params1.Frequency = [-50 50]; % MHz
Params1.tp = 0.200; % us
Params1.beta = 10;
Params1.Amplitude = PulseAmplitudes(1); % MHz

[t1_ideal,signal1_ideal] = pulse(Params1);
[offset1_ideal,M1_ideal] = exciteprofile(t1_ideal,signal1_ideal);

% Pulse 2
Params2 = Params1;
Params2.tp = 0.100; % us
Params2.Amplitude = PulseAmplitudes(2); % MHz
[t2_ideal,signal2_ideal] = pulse(Params2);
[offset2_ideal,M2_ideal] = exciteprofile(t2_ideal,signal2_ideal);

% Amplitude compression
% ---------------------------------------------------------------------------------------
Opt.Transmitter = 'simulate';

% Pulse 1
Params1.Amplitude = PulseAmplitudes(1)/max(PulseAmplitudes);
Sin =  0:0.001:1;
Gain = 85.5; % extrapolated output amplitude for max input amplitude in linear regime
Input1dBPoint = 0.60;
Sout = (Gain*Sin - (Gain*(1-10^(-1/20))/Input1dBPoint^2)*Sin.^3);
[~,ind] = min(abs(diff(Sout)));
Ssat = Sout(ind);
Sout(ind:end) = Ssat;
Params1.AmplifierResponse = [Sin; Sout]; % relative scale, MHz

[t1_compressed,signal1_compressed] = pulse(Params1,Opt);
[offset1_compressed,M1_compressed] = exciteprofile(t1_compressed,signal1_compressed);

% Pulse 2
Params2.Amplitude = PulseAmplitudes(2)/max(PulseAmplitudes);
Params2.AmplifierResponse = Params1.AmplifierResponse;

[t2_compressed,signal2_compressed] = pulse(Params2,Opt);
[offset2_compressed,M2_compressed] = exciteprofile(t2_compressed,signal2_compressed);

% Amplitude compensation
% ---------------------------------------------------------------------------------------
nu1_out = max(Params1.AmplifierResponse(2,:));
Opt.Transmitter = 'compensate';

% Pulse 1
Params1.Amplitude = nu1_out*Params1.Amplitude; % MHz
[t1_compensated,signal1_compensated] = pulse(Params1,Opt);

% Pulse 2
Params2.Amplitude = nu1_out*Params2.Amplitude; % MHz
[t2_compensated,signal2_compensated] = pulse(Params2,Opt);

% Effect of the transmitter on compensated pulses
% ---------------------------------------------------------------------------------------
Opt.Transmitter = 'simulate';

% Pulse 1
Params1_.tp = 0.200; % us
Params1_.I = real(signal1_compensated);
Params1_.Q = imag(signal1_compensated);

Params1_.Amplitude = max(abs(signal1_compensated));
Params1_.AmplifierResponse = Params1.AmplifierResponse;

[t1_corrected,signal1_corrected] = pulse(Params1_,Opt);
[offset1_corrected,M1_corrected] = exciteprofile(t1_corrected,signal1_corrected);

% Pulse 2
Params2_.tp = 0.100; % us
Params2_.I = real(signal2_compensated);
Params2_.Q = imag(signal2_compensated);

Params2_.Amplitude = max(abs(signal2_compensated));
Params2_.AmplifierResponse = Params2.AmplifierResponse;

[t2_corrected,signal2_corrected] = pulse(Params2_,Opt);
[offset2_corrected,M2_corrected] = exciteprofile(t2_corrected,signal2_corrected);

% Plot results
% ---------------------------------------------------------------------------------------
% Transmitter characterization
subplot(2,2,1)
hold on; box on; grid on;
title('Transmitter nonlinearity characterization')
plot(Params1.AmplifierResponse(1,:),Params1.AmplifierResponse(1,:)*Gain,'--k')
plot(Params1.AmplifierResponse(1,:),Params1.AmplifierResponse(2,:),'b')
line([1 1]*Input1dBPoint,ylim,'LineStyle','--','Color',[0.7 0.7 0.7])
line(xlim,[1 1]*Gain*Input1dBPoint*10^(-1/20),'LineStyle','--','Color',[0.7 0.7 0.7])
text(0.02,Gain*Input1dBPoint*10^(-1/20)+5,'1 dB compression','FontSize',7)
line(xlim,[1 1]*max(Params2_.AmplifierResponse(2,:)),'LineStyle','--','Color','k')
text(0.02,max(Params2_.AmplifierResponse(2,:))+5,'saturation','FontSize',7)
xlabel('Input amplitude (norm.)')
ylabel('Output amplitude (MHz)')

% Pulse shapes
subplot(2,2,2)
hold on; box on; grid on;
title('Pulses')
plot(t1_ideal,real(signal1_ideal),'Color',[0.7 0.7 0.7])
plot(t1_compressed,real(signal1_compressed),'r')
plot(t1_corrected,real(signal1_corrected),'b')
legend('ideal','compressed','corrected','Location','Best')
plot(t2_ideal+delta,real(signal2_ideal),'Color',[0.7 0.7 0.7])
plot(t2_compressed+delta,real(signal2_compressed),'r')
plot(t2_corrected+delta,real(signal2_corrected),'b')
axis([0 0.4 -65 65])
line(xlim,[0 0],'Color','k')
line(xlim,[1 1]*max(Params2_.AmplifierResponse(2,:)),'Color','k')
xlabel(['t (',char(181),'s)'])
ylabel('\nu_1 (MHz)')

% Compensated pulse shapes
subplot(2,2,3)
hold on; box on; grid on;
title('Compensated pulse shapes')
plot(t1_ideal,real(signal1_ideal)/max(real(signal2_ideal)),'Color',[0.7 0.7 0.7])
plot(t1_compensated,real(signal1_compensated),'b')
legend('ideal','compensated','Location','Best')
plot(t2_ideal+delta,real(signal2_ideal)/max(real(signal2_ideal)),'Color',[0.7 0.7 0.7])
plot(t2_compensated+delta,real(signal2_compensated),'b')
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
plot(offset1_corrected,M1_corrected(3,:),'b')
legend('ideal','compressed','corrected','Location','Best')
plot(offset2_ideal,M2_ideal(3,:),'--','Color',[0.7 0.7 0.7])
plot(offset2_compressed,M2_compressed(3,:),'--r')
plot(offset2_corrected,M2_corrected(3,:),'--b')
xlabel('\Delta\nu (MHz)')
ylabel('M_z/M_0')
