% Comparison of excitation profiles of a rectangular and sech/tanh pulse
% with the same maximum amplitude
%==========================================================================
% In this example, the excitation profiles of rectangular and sech/tanh pi 
% pulses with the same maximum amplitude are computed and the inversion
% profiles are compared.

clear, clf

% Sech/tanh pulse
%-------------------------------------------------------------
Par{1}.Type = 'sech/tanh'; % pulse shape
Par{1}.tp = 0.200; % pulse length, µs
Par{1}.Frequency = [-100 100]; % pulse frequency sweep range, MHz
Par{1}.beta = 8; % truncation parameter, used as (beta/tp)
Par{1}.Flip = pi; % pulse flip angle

Opt.Offsets = -200:1:200;

[t{1},IQ{1},modulation{1}] = pulse(Par{1},Opt);

% Rectangular pulse
%-------------------------------------------------------------
Par{2}.Type = 'rectangular';
% Use same amplitude as for the sech/tanh pi pulse
Par{2}.Amplitude = max(modulation{1}.A); % pulse amplitude, MHz
% Calculate the rectangular pulse length for a pi pulse
tp = pi/(2*pi*Par{2}.Amplitude); % pulse length = flip angle/amplitude, µs
Par{2}.tp = round(tp*1e3)/1e3;

[t{2},IQ{2}] = pulse(Par{2},Opt);

% Calculate excitation profiles
%-------------------------------------------------------------
offsets = -200:1:200; % offset axis in MHz
for ipulse = 1:2
  [offsets,M{ipulse}] = exciteprofile(t{ipulse},IQ{ipulse},offsets); 
end

% Plot
%--------------------------------------------------------------------
colI = [0 0 1];
colQ = [1 0 0];

subplot(3,1,1)
hold on; box on;
title([num2str(Par{2}.tp*1e3) ' ns rectangular pulse']);
plot(t{2},real(IQ{2}),'Color',colI)
plot(t{2},imag(IQ{2}),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
mA = max(modulation{1}.A);
axis([t{1}(1) t{1}(end) -1.1*mA 1.1*mA]);

subplot(3,1,2)
hold on; box on;
title([num2str(Par{1}.tp*1e3) ' ns sech/tanh pulse']);
plot(t{1},real(IQ{1}),'Color',colI)
plot(t{1},imag(IQ{1}),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
axis([t{1}(1) t{1}(end) -1.1*mA 1.1*mA]);

subplot(3,1,3)
hold on; box on;
title('Inversion profiles');
plot(offsets,M{2}(3,:),...
     offsets,M{1}(3,:))
xlabel('\Delta\nu (MHz)')
ylabel('M_z/M_0')
legend('rectangular','sech/tanh')
