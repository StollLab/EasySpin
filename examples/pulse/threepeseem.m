% manual simulation of three-pulse ESEEM (density matrix)
%==========================================================================
% Here we show how one can simulate a 3-pulse ESEEM spectrum
% manually, without using EasySpins' pulse simulation function
% saffron. This is an example of advanced usage.

clear, clf

n = 1000; % signal length
dT = .005; % step time [탎]
T0 = .300; % start T [탎]
tau = .100; % tau [탎]
FWHM = 20; % line width [MHz]

% spin Hamiltonian parameters.
A = 3.25;
B = 0.46;
wn = 6;

% spin operators
Sys = [1/2 1/2];
Sx = sop(Sys,'xe'); Sy = sop(Sys,'ye'); Sz = sop(Sys,'ze');
Ix = sop(Sys,'ex'); Iy = sop(Sys,'ey'); Iz = sop(Sys,'ez');

% The static Hamiltonian, on-resonance, rotating frame
Ham0 = A*Sz*Iz + B*Sz*Ix + wn*Iz;

offset = 2*FWHM*linspace(-1,1,500);
weights = gaussian(offset,0,FWHM);

% Propagation operators
Pulse90 = expm(-1i*pi/2*Sy);
Density = -Sz;

% Loop and sum over all offsets
signal = zeros(n,1);
for k=1:length(offset)
  Ham = offset(k)*Sz + Ham0;
  TauEvol = expm(-2i*pi*tau*Ham);
  T0Evol = expm(-2i*pi*T0*Ham);
  Prep = T0Evol*Pulse90*TauEvol*Pulse90;
  PrepDens = Prep*Density*Prep';
  Post = TauEvol*Pulse90;
  Detect = Post'*Sx*Post;
  signal = signal + weights(k)*evolve(PrepDens,Detect,Ham,n,dT,[1]);
end
signal = real(signal);
T = (0:n-1)*dT + T0;

subplot(2,1,1);
plot(T,signal,'b');
title('Three-pulse ESEEM time-domain signal');
xlabel('T [탎]');

subplot(2,1,2);
nn = 10*n;

alpha = 6;
KaiserWin = apowin('kai',length(signal),alpha).';

spec = fftshift(fft((signal-mean(signal)).*KaiserWin.',nn));
plot(fdaxis(dT,nn),abs(spec));
title('Three-pulse ESEEM magnitude spectrum');
xlim([0 1]*20);
xlabel('freq [MHz]');

