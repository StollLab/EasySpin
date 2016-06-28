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
Par{1}.BW = 180; % pulse bandwidth, MHz
Par{1}.beta = 8; % truncation parameter, used as (beta/tp)
Par{1}.Flip = pi; % pulse flip angle

Opt.Offsets = -200:1:200;

[t{1},y{1},p{1},m{1}] = pulse(Par{1},Opt);

% Rectangular pulse
%-------------------------------------------------------------
Par{2}.Type = 'rectangular';
% Use same amplitude as for the sech/tanh pi pulse
Par{2}.Amplitude = max(m{1}.A); % pulse amplitude, MHz
% Calculate the rectangular pulse length for a pi pulse
tp = pi/(2*pi*Par{2}.Amplitude); % pulse length = flip angle/amplitude, µs
Par{2}.tp = round(tp*1e3)/1e3;

[t{2},y{2},p{2}] = pulse(Par{2},Opt);

% Plot
%--------------------------------------------------------------------
colI = [0 0 1];
colQ = [1 0 0];

subplot(3,1,1)
hold on; box on;
title([num2str(Par{2}.tp*1e3) ' ns rectangular pulse']);
plot(t{2},real(y{2}),'Color',colI)
plot(t{2},imag(y{2}),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
mA = max(m{1}.A);
axis([t{1}(1) t{1}(end) -1.1*mA 1.1*mA]);

subplot(3,1,2)
hold on; box on;
title([num2str(Par{1}.tp*1e3) ' ns sech/tanh pulse']);
plot(t{1},real(y{1}),'Color',colI)
plot(t{1},imag(y{1}),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
axis([t{1}(1) t{1}(end) -1.1*mA 1.1*mA]);

subplot(3,1,3)
hold on; box on;
title('Inversion profiles');
plot(p{2}.offsets,p{2}.Mz,p{1}.offsets,p{1}.Mz)
xlabel('\Delta\nu (MHz)')
ylabel('M_z/M_0')
legend('rectangular','sech/tanh')
