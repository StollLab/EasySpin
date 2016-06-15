function PulsePropagator = sf_propagator(t,pulse,H0,Sx,Sy)

%
% Computes pulse propagator for arbitrarily shaped pulse
%
% Input: t     = pulse time axis
%        pulse = real and imaginary amplitude of the pulse
%        H0    = static Hamiltonian in correct frame
%        Sx    = spin operator Sx in correct frame
%        Sy    = spin operator Sy in correct frame
%

% time increment and pulse length
dt = t(2) - t(1);
tp = t(end) - t(1);

if min(pulse)==max(pulse) % rectangular pulse
  
  H = H0 + real(pulse(1))*Sx + imag(pulse(1))*Sy;
  PulsePropagator = expm(-2i*pi*H*tp);
  
else % amplitude/frequency modulated pulse
  
  for pt = 1:numel(pulse)-1
    
    H = H0 + real(pulse(pt))*Sx + imag(pulse(pt))*Sy;
    if pt==1
      PulsePropagator = expm(-2i*pi*H*dt);
    else
      PulsePropagator = expm(-2i*pi*H*dt)*PulsePropagator;
    end
    
  end
  
end
