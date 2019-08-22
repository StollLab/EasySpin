% Resonator bandwidth compensation by chirp rate adaptation
%==========================================================================
% In this example, the frequency and amplitude modulation functions of a 
% sech/tanh pulse are adapted to account for the limited resonator band-
% width and achieve offset-independent adiabaticity.
%
% This method is described in Doll, A., Pribitzer, S., Tschaggelar, R.,
% Jeschke, G., J. Magn. Reson. 230, 27-39 (2013)
% https://doi/org/10.1016/j.jmr.2013.01.002

clear, clf, clc

% Define sech/tanh pulse
%-------------------------------------------------------------------------------
Par.Type = 'sech/tanh'; % pulse shape
Par.tp = 0.200; % pulse length, microseconds
Par.Frequency = [-80 80]; % pulse frequency sweep range, MHz
Par.beta = 8; % truncation parameter, used as (beta/tp)
Par.Flip = pi; % pulse flip angle

% Calculate pulse without resonator compensation
%-------------------------------------------------------------------------------
[t,IQ,modulation] = pulse(Par);

% Recalculate pulse with resonator compensation (amplitude and frequency)
%-------------------------------------------------------------------------------
Par.mwFreq = 33.85;  % microwave frequency, GHz
Par.ResonatorFrequency = 33.9;  % resonator center frequency, GHz
Par.ResonatorQL = 200;     % loaded Q-factor of resonator
[t_comp,IQ_comp,modulation_comp] = pulse(Par);

% Plotting
%-------------------------------------------------------------------------------
colI = [0 0 1]; colIfaded = 0.25*colI + 0.75*[1 1 1];
colQ = [1 0 0]; colQfaded = 0.25*colQ + 0.75*[1 1 1];

% Transfer function (for plotting)
f = Par.ResonatorFrequency + linspace(-1,1,1001)*0.5; % GHz
H = 1./(1+1i*Par.ResonatorQL*(f/Par.ResonatorFrequency-Par.ResonatorFrequency./f));

maxAmplitude = max(modulation.A);

subplot(2,2,1)
hold on; box on;
title('resonator transfer function')
pulseFreqRange = Par.mwFreq+Par.Frequency*1e-3;
patch(pulseFreqRange([1 1 2 2]),[0 1 1 0]*1.5,[1 1 1]*0.85,'EdgeColor','none');
text(pulseFreqRange(1),0.5,'pulse ','FontSize',8,'Color',[1 1 1]*0.5,'HorizontalA','right');
line([1 1]*Par.mwFreq,ylim,'Color','k')
text(Par.mwFreq+0.01,1.25,'mwFreq','FontSize',8)
plot(f,real(H),'b')
axis tight
set(gca,'Layer','top')
xlabel('frequency (GHz)')
ylabel('transfer function')

subplot(2,2,3)
hold on; box on;
title('frequency modulation');
plot(t,modulation.freq,'Color',colIfaded)
plot(t_comp,modulation_comp.freq,'Color',colI)
xlabel('t (\mus)')
ylabel('frequency (MHz)')
axis tight
ylim(1.1*Par.Frequency);

subplot(2,2,4)
hold on; box on;
title('ampltiude modulation');
plot(t,modulation.A,'Color',colIfaded)
plot(t_comp,modulation_comp.A,'Color',colI)
xlabel('t (\mus)')
ylabel('amplitude (MHz)')
axis tight
ylim([0 1]*1.1*maxAmplitude);

subplot(2,2,2)
hold on; box on;
title('sech/tanh pulse w/o and w/ resonator BW compensation');
plot(t,real(IQ),'Color',colIfaded)
plot(t,imag(IQ),'Color',colQfaded)
plot(t_comp,real(IQ_comp),'Color',colI)
plot(t_comp,imag(IQ_comp),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
axis tight
ylim([-1 1]*1.1*maxAmplitude);
