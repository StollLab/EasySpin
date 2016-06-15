% manual two-pulse ESEEM simulation (density matrix)
%==========================================================================
% Here we show how one can manually simulate a 2-pulse ESEEM
% spectrum using the density matrix formalism. This is an
% example of advanced usage.

clear, clf

n = 1000; % signal length
dtau = .005; % step time [µs]
FWHM = 20; % line width [MHz]

% spin Hamiltonian parameters.
A = 3.25;
B = 0.46;
wn = 6;

% spin operators
Sys = [1/2 1/2];
Sx = sop(Sys,'xe'); Sy = sop(Sys,'ye'); Sz = sop(Sys,'ze');
Ix = sop(Sys,'ex'); Iy = sop(Sys,'ey'); Iz = sop(Sys,'ez');

% The static Hamiltonian
Ham0 = A*Sz*Iz + B*Sz*Ix + wn*Iz;

offset = 1.2*FWHM*linspace(-1,1,200);
weights = gaussian(offset,0,FWHM);

% Propagation operators
Pulse = expm(-1i*pi/2*Sx);
Dens = Pulse*(Sz)*Pulse';
Mix = expm(-1i*pi*Sx);


% Loop and sum over all offsets
s = zeros(n,1);
for k=1:length(offset)
  Ham = offset(k)*Sz + Ham0;
  s = s + weights(k)*real(evolve(Dens,Sy,Ham,n,dtau,[1 1],{Mix}));
end

subplot(2,1,1);
plot((0:n-1)*dtau,s,'b');
title('Two-pulse ESEEM time-domain signal');
xlabel('tau [µs]');

subplot(2,1,2);
nn = 10*n;

alpha = 6;
KaiserWin = apowin('kai',length(s),alpha).';

spec = fftshift(fft((s-mean(s)).*KaiserWin.',nn));
plot(fdaxis(dtau,nn),abs(spec));
title('Two-pulse ESEEM magnitude spectrum');
xlim([0 1]*20);
xlabel('freq [MHz]');
