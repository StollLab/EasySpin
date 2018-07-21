% Simulation of a two-pulse echo using the density matrix
%==========================================================================
% Here we show how one can simulate a 2-pulse echo transient
% using the density matrix formalism. This is an example of
% advanced usage.

clear, clf

n = 300; % signal length
dt = 0.005; % step time [탎]
FWHM = 20; % line width [MHz]
tau = 0.7; % interpulse distance [탎]

Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

offset = 1.2*FWHM*linspace(-1,1,200);
weights = gaussian(offset,0,FWHM);
tp = 0:0.020:0.180; % pulse length [탎]

% We compute a set of 2-pulse-echo signals, 
% differing in the length of the pulses used.
% The static Hamiltonian is not neglected during
% the pulse.

flipAngle = pi/2;
signal = zeros(n,length(tp));
for p = 1:length(tp)
  for k = 1:length(offset)
    Ham0 = offset(k)*Sz;
    Pulse = expm(-1i*(flipAngle*Sx+2*pi*tp(p)*Ham0));
    TauEvol = expm(-2i*pi*tau*Ham0);
    U = Pulse^2*TauEvol*Pulse;
    signal(:,p) = signal(:,p) + ...
      weights(k)*real(evolve(U*Sz*U',Sy,Ham0,n,dt));
  end
end

% Result: The longer the pulses, the broader the
% echo.

plot((0:n-1)*dt,signal,'b');
title('The shape of the primary echo depending on the pulse length');
xlabel('time after the pi pulse [탎]');
ylabel('echo signal');

% end
