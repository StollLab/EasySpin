% Time-domain Mims ENDOR of a 1H (saffron)
%==========================================================================
% Powder time-domain ENDOR spectrum of an S=1/2 with a weakly coupled proton.
% Instead of sweeping the frequency of an RF pi pulse, time-domain ENDOR
% splits it into two RF pi/2 pulses at a fixed frequency and increments the
% delay t between them:
%   pi/2 - tau - pi/2 - T - RF pi/2 - t - RF pi/2 - T - pi/2 - tau - echo
% The first RF pulse creates nuclear coherence, which evolves during t, and
% the second RF pulse converts it back into nuclear populations, which the
% stimulated echo detects. The echo oscillates with the differences between
% the nuclear transition frequencies and the RF frequency, and its Fourier
% transform gives the ENDOR spectrum. Unlike frequency-domain ENDOR, the
% line widths are not broadened by the RF pulse length.
%
% Here the RF frequency is set to the 1H Larmor frequency, and the
% spectrum is the Pake pattern of the electron-proton dipolar coupling,
% centered at the Larmor frequency. The echo is a real signal, so it gives
% only the distance of each line from the RF frequency, not on which side
% the line is. Its Fourier transform is therefore symmetric about the RF
% frequency. This is correct here, since the Pake pattern is symmetric
% about the Larmor frequency. With the RF frequency elsewhere, each line
% would also appear mirrored on the other side.
%
% As in frequency-domain Mims ENDOR, the intensity is zero for a coupling
% of zero (blind spot). This causes the dip at the center of the pattern,
% from orientations near the magic angle, where the dipolar coupling
% vanishes.
%
% The RF pulses must be strong enough to excite the whole ENDOR spectrum.
% saffron simulates this experiment with the full method, which takes
% some time.
%
% References:
%   P. Höfer, A. Grupp, M. Mehring, High-resolution time-domain
%     electron-nuclear-sublevel spectroscopy by pulsed coherence transfer,
%     Phys. Rev. A 33, 3519-3522 (1986), doi:10.1103/PhysRevA.33.3519
%   H. Cho, Stimulated echo time domain electron nuclear double resonance,
%     J. Chem. Phys. 94, 2482-2495 (1991), doi:10.1063/1.459872
%   A. Kehl, L. Sielaff, L. Remmel, M. L. Rämisch, M. Bennati, A. Meyer,
%     Frequency and time domain 19F ENDOR spectroscopy: role of nuclear
%     dipolar couplings to determine distance distributions,
%     Phys. Chem. Chem. Phys. 27, 1415-1425 (2025), doi:10.1039/d4cp04443f

clear, clc

% Spin system: S=1/2 with one 1H
Sys.g = gfree;
Sys.Nucs = '1H';
Tdip = 0.2;               % dipolar coupling (ca 7.3Å distance), MHz
Sys.A = [-Tdip 2*Tdip];   % axial hyperfine tensor, MHz

% Pulses
p90.Flip = pi/2;          % ideal microwave pulse flip angle, rad
pRF.Channel = 'rf';       % RF pulse
pRF.Flip = pi/2;          % RF pulse flip angle, rad
pRF.tp = 1;               % RF pulse length, µs

% Pulse sequence
tau = 1.2;                % µs
T = 1;                    % µs
t0 = 0;                   % µs
Exp.Sequence = {p90 tau p90 T pRF t0 pRF T p90 tau};

Exp.mwFreq = 9.8;         % GHz
Exp.Field = mhz2mt(Exp.mwFreq*1e3);      % mT, on resonance for g = 2.0023
Exp.rfFreq = larmorfrq('1H',Exp.Field);  % RF frequency = 1H Larmor frequency, MHz

% Electron coherence transfer pathway and detection phase
Exp.Pathways = '+000-';   % pathway for 3-pulse echo
Exp.DetPhase = -pi/2;     % detection phase for real positive echo, rad

% Increment the delay between the two RF pulses (third delay)
Exp.Dim1 = {'d3', 0.5};   % time increment, µs
Exp.nPoints = 101;

Opt.GridSize = 31;        % powder integration grid

% Simulation
[t,V,info] = saffron(Sys,Exp,Opt);
nuRF = Exp.rfFreq + info.f;   % RF frequency axis, MHz
spc = abs(info.fd);           % ENDOR spectrum (magnitude of Fourier transform)

% Plot time trace
subplot(2,1,1)
plot(t,real(V));
axis tight
xlabel('t (µs)');
ylabel('echo amplitude');
title('Time-domain Mims ENDOR');

% Plot ENDOR spectrum
subplot(2,1,2)
plot(nuRF,spc);
axis tight
xline(Exp.rfFreq,'--');
xlabel('RF frequency (MHz)');
ylabel('ENDOR intensity');
title('Fourier transform');
