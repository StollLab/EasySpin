% cw EPR as a multiphoton experiment
%==========================================================================
% This example illustrates the fact that the apparently simple cw EPR
% spectrum with a single line in reality is the result of overlap of many
% modulation sidebands.
% See Kaelin et al, J.Magn.Reson. 160, 166-182 (2003)

clear, clf

% absorption spectrum, first harmonic, in-phase with modulation
mwPhase = 0;  % mw phase (0 absorption, pi/2 dispersion)
Harmonic = 1; % detection harmonic (0, 1, 2, etc)
modPhase = 0; % modulation phase (0 in-phase, pi/2 out-of-phase)

fmod = 1; % modulation frequency (MHz)
Amod = 1;  % peak-to-peak modulation amplitude (mT)
lwpp = 0.1; % peak-to-peak line width (mT)

wrf = 2*pi*fmod; % angular modulation frequency (2*pi*MHz)
w2 = 2*pi*mt2mhz(Amod)/2; % Amplitude of the modulation field (2*pi*MHz)
T2 = 1/mt2mhz(lwpp); % Relaxation time (us)

nPoints = 10000;
fieldmax = max(Amod*2,lwpp*8);
field = fieldmax*linspace(-1,1,nPoints);
w = 2*pi*mt2mhz(field);

% The following implements Eqs.(59) and (60) from the Kaelin paper,
% disregaring the saturation term in the denominators.
% In the paper, Eq.(59) is mislabeled as Eq.(62). Also, in that
% equation, factors j are missing in the cos and sin functions.
spec = 0;
j = Harmonic;
beta = w2/wrf; % modulation parameter
kmax = ceil(max(1.1*beta,beta+30));
for k = -kmax:kmax
  jm = besselj(-k,beta);
  j1 = besselj(-k+j,beta);
  j2 = besselj(-k-j,beta);
  Acos = jm*(j1+j2); % Eqs.(59) and (60), cos(j*wrf*t) term
  Asin = jm*(j1-j2); % Eqs.(59) and (60), sin(j*wrf*t) term
  denom = 1 + (w-k*wrf).^2*T2^2; % disregard saturation term
  shape1 = T2./denom;
  shape2 = (w-k*wrf).*T2^2./denom;
  Abs0  = +shape1*Acos;
  Abs90 = -shape2*Asin;
  Dis0  = -shape2*Acos;
  Dis90 = -shape1*Asin;
  spec = spec + ...
         (Abs0*cos(modPhase) + Abs90*sin(modPhase))*cos(mwPhase) + ...
         (Dis0*cos(modPhase) + Dis90*sin(modPhase))*sin(mwPhase);
end

plot(field,spec);
xlabel('field offset (mT)');
ylabel('amplitude (arb.u.)');
