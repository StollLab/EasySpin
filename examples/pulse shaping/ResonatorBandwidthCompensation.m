% Resonator bandwidth compensation by chirp rate adaptation (as described 
% in Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G., J. Magn. Reson.
% 230, 27-39 (2013), DOI: 10.1016/j.jmr.2013.01.002)
%==========================================================================
% In this example, the frequency and amplitude modulation functions of a 
% sech/tanh pulse are adapted to account for the limited resonator band-
% width and achieve offset-independent adiabaticity.

clear, clf

% Sech/tanh pulse
%--------------------------------------------------------------------
Par.Type = 'sech/tanh'; % pulse shape
Par.tp = 0.200; % pulse length, µs
Par.Frequency = [-80 80]; % pulse frequency sweep range, MHz
Par.beta = 8; % truncation parameter, used as (beta/tp)
Par.Flip = pi; % pulse flip angle

[t{1},IQ{1},modulation{1}] = pulse(Par);

% Resonator properties
%--------------------------------------------------------------------
Par.mwFreq = 33.85;  % microwave frequency, GHz
Par.ResonatorFrequency = 33.9;  % resonator center frequency, GHz
Par.ResonatorQL = 200;     % loaded Q-factor of resonator

% Recalculate pulse with adapted sweep rate 
%--------------------------------------------------------------------
[t{2},IQ{2},modulation{2}] = pulse(Par);

% Plot
%--------------------------------------------------------------------
colI = [0 0 1];
colQ = [1 0 0];

% Transfer function (for plot)
f = 33.5:0.001:34.3;
H = 1./(1+1i*Par.ResonatorQL*(f/Par.ResonatorFrequency-Par.ResonatorFrequency./f));

subplot(3,1,1)
hold on; box on;
title('Resonator transfer function')
patch([Par.mwFreq+Par.Frequency(1)*1e-3 Par.mwFreq+Par.Frequency(1)*1e-3 ...
       Par.mwFreq+Par.Frequency(2)*1e-3 Par.mwFreq+Par.Frequency(2)*1e-3],...
      [0 1.5 1.5 0],[1 1 1]*0.85,'EdgeColor','none');
text(Par.mwFreq+Par.Frequency(1)*1e-3-0.065,0.5,'pulse','FontSize',8,'Color',[1 1 1]*0.85)
line([1 1]*Par.mwFreq,ylim,'Color','k')
text(Par.mwFreq+0.01,1.25,'mwFreq','FontSize',8)
plot(f,real(H),'b')
axis tight
set(gca,'Layer','top')
xlabel('Frequency (GHz)')
ylabel('Transfer function')

subplot(3,1,2)
hold on; box on;
title('frequency modulation w/o and w/ resonator BW compensation');
plot(t{2},modulation{1}.freq,'Color',0.25*colI+0.75*ones(1,3))
plot(t{2},modulation{2}.freq,'Color',colI)
xlabel('t (\mus)')
ylabel('frequency (MHz)')
axis([t{1}(1) t{1}(end) 1.1*Par.Frequency]);

subplot(3,1,3)
hold on; box on;
title('sech/tanh pulse w/o and w/ resonator BW compensation');
plot(t{1},real(IQ{1}),'Color',0.25*colI+0.75*ones(1,3))
plot(t{1},imag(IQ{1}),'Color',0.25*colQ+0.75*ones(1,3))
plot(t{2},real(IQ{2}),'Color',colI)
plot(t{2},imag(IQ{2}),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
mA = max(modulation{1}.A);
axis([t{1}(1) t{1}(end) -1.1*mA 1.1*mA]);
