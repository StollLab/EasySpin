% Illustration of the effect of the resonator on pulse shapes and pulse
% compensation based on a transfer function
%==========================================================================
% In this example, the effect of the resonator on a rectangular pulse is
% investigated and the transfer function is used to modify the pulse shape
% to compensate for the distortions in the resonator. By simulating the 
% effect of the resonator on the compensated pulse shape the success of the
% correction is shown.

clear, clc, clf

% Resonator properties
% ---------------------------------------------------------------------------------------
ResonatorFrequency = 9.50; % GHz
ResonatorQL = 300;
mwFreq = 9.49; % GHz
% Transfer function (for plot)
f = 9.2:0.001:9.8;
H = 1./(1+1i*ResonatorQL*(f/ResonatorFrequency-ResonatorFrequency./f));

% Ideal pulses
% ---------------------------------------------------------------------------------------
Params1.Type = 'quartersin';
Params1.trise = 0.004; % ns
Params1.Frequency = 0; % MHz
Params1.tp = 0.100; % us
Params1.Amplitude = 15; % MHz

[t1_ideal,signal1_ideal] = pulse(Params1);

% Pulse 2
Params2 = Params1;
Params2.Frequency = -20; % MHz
[t2_ideal,signal2_ideal] = pulse(Params2);

% Pulse 3
Params3.Type = 'sech/tanh';
Params3.tp = 0.200; % us
Params3.beta = 10;
Params3.Frequency = [-50 50]; % MHz
Params3.Amplitude = 15; % MHz
[t3_ideal,signal3_ideal] = pulse(Params3);


% Effect of the resonator
% ---------------------------------------------------------------------------------------
% Pulse 1
[t1_dist,signal1_dist] = resonator(t1_ideal,signal1_ideal,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Pulse 2
[t2_dist,signal2_dist] = resonator(t2_ideal,signal2_ideal,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Pulse 3
[t3_dist,signal3_dist] = resonator(t3_ideal,signal3_ideal,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Resonator compensation
% ---------------------------------------------------------------------------------------
Opt.Resonator = 'compensate';

% Pulse 1
[t1_compensated,signal1_compensated] = resonator(t1_ideal,signal1_ideal,mwFreq,ResonatorFrequency,ResonatorQL,'compensate');

% Pulse 2
[t2_compensated,signal2_compensated] = resonator(t2_ideal,signal2_ideal,mwFreq,ResonatorFrequency,ResonatorQL,'compensate');

% Pulse 3
[t3_compensated,signal3_compensated] = resonator(t3_ideal,signal3_ideal,mwFreq,ResonatorFrequency,ResonatorQL,'compensate');

% Effect of the resonator on compensated pulses
% ---------------------------------------------------------------------------------------

% Pulse 1
[t1_corrected,signal1_corrected] = resonator(t1_compensated,signal1_compensated,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Pulse 2
[t2_corrected,signal2_corrected] = resonator(t2_compensated,signal2_compensated,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Pulse 3
[t3_corrected,signal3_corrected] = resonator(t3_compensated,signal3_compensated,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');


% Plot results
% ---------------------------------------------------------------------------------------
% Resonator characterization
subplot(3,1,1)
hold on; box on;
title('Resonator transfer function')
patch([mwFreq+Params3.Frequency(1)*1e-3 mwFreq+Params3.Frequency(1)*1e-3 ...
       mwFreq+Params3.Frequency(2)*1e-3 mwFreq+Params3.Frequency(2)*1e-3],...
      [-1 1 1 -1],[1 1 1]*0.85,'EdgeColor','none');
text(mwFreq+Params3.Frequency(1)*1e-3-0.065,-0.4,'pulse 3','FontSize',8,'Color',[1 1 1]*0.85)
line([1 1]*mwFreq,ylim,'Color','k')
text(mwFreq+0.01,-0.75,'mwFreq','FontSize',8)
line([1 1]*mwFreq+Params2.Frequency*1e-3,ylim,'Color',[1 1 1]*0.6)
text(mwFreq+Params2.Frequency*1e-3-0.065,-0.65,'pulse 2','FontSize',8,'Color',[1 1 1]*0.6)
plot(f,real(H),'b')
plot(f,imag(H),'r')
axis tight
set(gca,'Layer','top')
xlabel('Frequency (GHz)')
ylabel('Transfer function')

% Pulse shapes
subplot(3,1,2:3)
delta = 0.200;
hold on; box on;
title('Pulse shapes')
line([0 2*Params1.tp+2*delta+0.050],[0 0],'Color','k')
line([0 2*Params1.tp+2*delta+0.050],[0 0]-3*Params1.Amplitude,'Color','k')
line([0 2*Params1.tp+2*delta+0.050],[0 0]-7*Params1.Amplitude,'Color','k')
line([0 2*Params1.tp+2*delta+0.050],[0 0]-11*Params1.Amplitude,'Color','k')
% Ideal pulse shapes
text(0.005,1.5*Params1.Amplitude,'Ideal pulses','FontSize',8)
plot(t1_ideal,real(signal1_ideal),'Color','b')
plot(t1_ideal,imag(signal1_ideal),'Color','r')
plot(t2_ideal+delta,real(signal2_ideal),'Color','b')
plot(t2_ideal+delta,imag(signal2_ideal),'Color','r')
plot(t3_ideal+2*delta,real(signal3_ideal),'Color','b')
plot(t3_ideal+2*delta,imag(signal3_ideal),'Color','r')

% Pulses distorted by the resonator
text(0.005,-1.5*Params1.Amplitude,'Pulses distorted by the resonator','FontSize',8)
plot(t1_dist,real(signal1_dist)-3*Params1.Amplitude,'Color','b')
plot(t1_dist,imag(signal1_dist)-3*Params1.Amplitude,'Color','r')
plot(t2_dist+delta,real(signal2_dist)-3*Params1.Amplitude,'Color','b')
plot(t2_dist+delta,imag(signal2_dist)-3*Params1.Amplitude,'Color','r')
plot(t3_dist+2*delta,real(signal3_dist)-3*Params1.Amplitude,'Color','b')
plot(t3_dist+2*delta,imag(signal3_dist)-3*Params1.Amplitude,'Color','r')

% Compensated pulse shapes
text(0.005,-4.5*Params1.Amplitude,'Compensated pulse shapes','FontSize',8)
plot(t1_compensated,real(signal1_compensated)-7*Params1.Amplitude,'Color','b')
plot(t1_compensated,imag(signal1_compensated)-7*Params1.Amplitude,'Color','r')
plot(t2_compensated+delta,real(signal2_compensated)-7*Params1.Amplitude,'Color','b')
plot(t2_compensated+delta,imag(signal2_compensated)-7*Params1.Amplitude,'Color','r')
plot(t3_compensated+2*delta,real(signal3_compensated)-7*Params1.Amplitude,'Color','b')
plot(t3_compensated+2*delta,imag(signal3_compensated)-7*Params1.Amplitude,'Color','r')

% Test
text(0.005,-9.5*Params1.Amplitude,'Effect of the resonator on the compensated pulses','FontSize',8)
plot(t1_corrected,real(signal1_corrected)-11*Params1.Amplitude,'Color','b')
plot(t1_corrected,imag(signal1_corrected)-11*Params1.Amplitude,'Color','r')
plot(t2_corrected+delta,real(signal2_corrected)-11*Params1.Amplitude,'Color','b')
plot(t2_corrected+delta,imag(signal2_corrected)-11*Params1.Amplitude,'Color','r')
plot(t3_corrected+2*delta,real(signal3_corrected)-11*Params1.Amplitude,'Color','b')
plot(t3_corrected+2*delta,imag(signal3_corrected)-11*Params1.Amplitude,'Color','r')

xlim([0 2*Params1.tp+2*delta+0.050])
xlabel(['t (',char(181),'s)'])
ylabel('\nu_1 (MHz)')
ylim([-13 2]*Params1.Amplitude)
set(gca,'YTick',0:10:20)
