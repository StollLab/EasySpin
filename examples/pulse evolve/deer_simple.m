% DEER from a simple system with two spins-1/2 (powder)
%======================================================================
% This example uses density matrix dynamics to simulate the DEER
% time trace for a simple two-spin-1/2 system.

clear, clc, clf

% Parameters
%----------------------------------------------------------------------
% Offset frequencies at given magnetic field
Nu1 = -130; % offset frequency probe spin, MHz
Nu2 = +30; % offset frequency pump spin, MHz
g1 = 2; % g-factor probe spin (only used to calculate dipole coupling)
g2 = 2; % g-factor pump spin (only used to calculate dipole coupling)

r = 3.5; % distance, nm

tau1 = 0.100;  % us
tau2 = 4;      % us
nPoints = 100; % number of points in the time domain
nAngles = 200; % number of orientations between theta = 0 and pi/2

% Preparations
%----------------------------------------------------------------------

t = linspace(0,tau2,nPoints); % us
V = zeros(nPoints,1);

theta = linspace(0,pi/2,nAngles);
weight = sin(theta);
weight = weight/sum(weight);
nuperp = (mu0/4/pi)*bmagn^2*g1*g2/(r*1e-9)^3/planck/1e6; % MHz
D = nuperp*(3*cos(theta).^2-1); % MHz

% Spin operators: 1 = probe spin, 2 = pump spin
[S1x,S1y,S1z] = sop([1/2 1/2],'x1','y1','z1');
[S2x,S2y,S2z] = sop([1/2 1/2],'x2','y2','z2');

% Pulse propagators
Probe90 = expm(-1i*pi/2*S1x); % 90-degree pulse on probe spins
Probe180 = Probe90^2;         % 180-degree pulse on probe spins
Pump180 = expm(-1i*pi*S2x);   % 180-degree pulse on pump spins

sigma_eq = -S1z-S2z;  % thermal equilibrium density
DetectOp = S1y; % detection operator

% Powder average
%----------------------------------------------------------------------
for itheta = 1:numel(theta)

  % Rotating-frame Hamiltonian for given orientation
  H = Nu1*S1z + Nu2*S2z + D(itheta)*S1z*S2z; % MHz
  
  % Propagate density from start to time of two-pulse echo maximum
  Utau = expm(-2i*pi*H*tau1);
  Prep = Utau*Probe180*Utau*Probe90;
  sigma_prep = Prep*sigma_eq*Prep';
  
  for it = 1:numel(t)
    % Calculate propagator for the rest of the sequence
    U1 = expm(-2i*pi*H*t(it)); % pre-pump
    U2 = expm(-2i*pi*H*(tau2-t(it))); % post-pump
    Prop = U2*U1*Probe180*U2*Pump180*U1; % total propagator
    
    % Propagate density
    sigma = Prop*sigma_prep*Prop';
    
    % Calculate expectation value and accumulate
    V(it) = V(it) + weight(itheta)*trace(DetectOp*sigma);
  end
  
end
V = real(V); % remove numerical noise in imaginary part (should be zero)

plot(t,V);
