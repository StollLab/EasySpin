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
Par{1}.BW = 200; % pulse bandwidth, MHz
Par{1}.nwurst = 8; % truncation parameter, used as (beta/tp)
Par{1}.Flip = pi/2; % pulse flip angle

[t{1},y{1},p{1}] = pulse(Par{1});

% pi pulse
%-------------------------------------------------------------
Par{2} = Par{1}; % same pulse shape parameters as the pi/2 pulse
Par{2}.tp = Par{1}.tp/2; % pulse length, µs
Par{2}.Flip = pi; % pulse flip angle

[t{2},y{2},p{2}] = pulse(Par{2});

% Plot
%--------------------------------------------------------------------
colI = [0 0 1];
colQ = [1 0 0];

subplot(2,2,1)
hold on; box on;
title('\pi/2 pulse');
plot(t{1},real(y{1}),'Color',colI)
plot(t{1},imag(y{1}),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
mA = max(abs(([y{1} y{2}])));
axis([t{1}(1) t{1}(end) -1.1*mA 1.1*mA]);

subplot(2,2,2)
hold on; box on;
title('\pi pulse');
plot(t{2},real(y{2}),'Color',colI)
plot(t{2},imag(y{2}),'Color',colQ)
xlabel('t (\mus)')
ylabel('\nu_1 (MHz)')
axis([t{1}(1) t{1}(end) -1.1*mA 1.1*mA]);

subplot(2,2,3)
hold on; box on;
title('Inversion profiles');
plot(p{1}.offsets,p{1}.Mx,p{1}.offsets,p{1}.My,p{1}.offsets,p{1}.Mz)
xlabel('\Delta\nu (MHz)')
ylabel('M_i/M_0')
legend('x','y','z')
axis tight

subplot(2,2,4)
hold on; box on;
title('Inversion profiles');
plot(p{2}.offsets,p{2}.Mx,p{2}.offsets,p{2}.My,p{2}.offsets,p{2}.Mz)
xlabel('\Delta\nu (MHz)')
ylabel('M_i/M_0')
legend('x','y','z')
axis tight