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
  
  if all(size(H)==2)
    % Fast matrix exponential for a 2x2 matrix
    M = -2i*pi*H;
    s = 0.5*(M(1,1)+M(2,2));
    q = sqrt(-det(M-s*eye(2)));
    if abs(q)<1e-10
      PulsePropagator = exp(s*tp)*((1-s)*eye(2) + M);
    else
      PulsePropagator = exp(s*tp)*((cosh(q*tp)-s*(sinh(q*tp)/q))*eye(2) + (sinh(q*tp)/q)*M);
    end
  else
    PulsePropagator = expm(-2i*pi*tp*H);
  end
  
else % amplitude/frequency modulated pulse
  
  PulsePropagator = eye(size(H0));
  for pt = 1:numel(pulse)-1
    
    H = H0 + real(pulse(pt))*Sx + imag(pulse(pt))*Sy;
    
    if all(size(H)==2)
      % Fast matrix exponential for a 2x2 matrix
      M = -2i*pi*H;
      s = 0.5*(M(1,1)+M(2,2));
      q = sqrt(-det(M-s*eye(2)));
      if abs(q)<1e-10
        U = exp(s*dt)*((1-s)*eye(2) + M);
      else
        U = exp(s*dt)*((cosh(q*dt)-s*(sinh(q*dt)/q))*eye(2) + (sinh(q*dt)/q)*M);
      end
      PulsePropagator = U*PulsePropagator;
    else
      PulsePropagator = expm(-2i*pi*dt*H)*PulsePropagator;
    end
    
  end
  
end
