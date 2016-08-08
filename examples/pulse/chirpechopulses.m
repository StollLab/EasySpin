% WURST chirp pi/2 and pi pulses
%==========================================================================
% In this example, WURST chirp pulses with flip angles of pi/2 and pi are
% computed and their excitation profiles are compared. The length of the pi
% pulse corresponds to half of the length of the pi/2 pulse following the
% Böhlen-Bodenhausen scheme for compensation of the nonlinear phase
% dispersion. The resulting pulses can be used to simulate a broadband
% two-pulse echo.

clear, clf

% pi/2 pulse
%-------------------------------------------------------------
Par{1}.Type = 'WURST/linear'; % pulse shape, linear chirp with WURST amplitude modulation
Par{1}.tp = 0.200; % pulse length, µs
Par{1}.Frequency = [-100 100]; % pulse bandwidth, MHz
Par{1}.nwurst = 8; % truncation parameter, used as (beta/tp)
Par{1}.Flip = pi/2; % pulse flip angle

[t{1},IQ{1}] = pulse(Par{1});

% pi pulse
%-------------------------------------------------------------
Par{2} = Par{1}; % same pulse shape parameters as the pi/2 pulse
Par{2}.tp = Par{1}.tp/2; % pulse length, µs
Par{2}.Flip = pi; % pulse flip angle

[t{2},IQ{2}] = pulse(Par{2});

% Calculate excitation profiles
%-------------------------------------------------------------
offsets = -150:1:150; % offset axis in MHz
for ipulse = 1:2
  
  [offsets,M{ipulse}] = exciteprofile(t{ipulse},IQ{ipulse},offsets);
  
end

% Plot
%--------------------------------------------------------------------
colI = [0 0 1];
colQ = [1 0 0];

subplot(2,2,1)
hold on; box on;
title('\pi/2 pulse');
plot(t{1},real(IQ{1}),'Color',colI)
plot(t{1},imag(IQ{1}),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
mA = max(abs(([IQ{1} IQ{2}])));
axis([t{1}(1) t{1}(end) -1.1*mA 1.1*mA]);

subplot(2,2,2)
hold on; box on;
title('\pi pulse');
plot(t{2},real(IQ{2}),'Color',colI)
plot(t{2},imag(IQ{2}),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
axis([t{1}(1) t{1}(end) -1.1*mA 1.1*mA]);

subplot(2,2,3)
hold on; box on;
title('Inversion profiles');
plot(offsets,M{1}(1,:),...
     offsets,M{1}(2,:),...
     offsets,M{1}(3,:))
xlabel('\Delta\nu (MHz)')
ylabel('M_i/M_0')
legend('x','y','z')
axis tight

subplot(2,2,4)
hold on; box on;
title('Inversion profiles');
plot(offsets,M{2}(1,:),...
     offsets,M{2}(2,:),...
     offsets,M{2}(3,:))
xlabel('\Delta\nu (MHz)')
ylabel('M_i/M_0')
legend('x','y','z')
axis tight