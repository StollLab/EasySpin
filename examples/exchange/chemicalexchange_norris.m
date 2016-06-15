% Multi-site chemical exchange
%===========================================================================
% Implements the Norris formula for multi-site exchange.
% The average lifetime in each site is tau, and the system
% can jump from each site to each other site.

% J.R.Norris, Chem.Phys.Lett. 1 (1967) 333-334
% Rapid Computation of Magnetic Resonance Lineshapes for Exchange Among Many Sites 

clear, clc

om0 = [1 5 12 19]; % resonance frequencies
p = [1 1 2 1.5]; % populations (Norris denotes them with n)
T2 = 10*[1 1 1 1]; % transverse relaxation time, potentially site-dependent
tau = 0.5; % average lifetime

om = linspace(-10,30,10001); % frequency axis
p = p/sum(p);

Fsum = 0;
N = numel(om0); % number of sites
for k = 1:N
  Fsum = Fsum + p(k)./(1i*(om-om0(k))-1/T2(k)-1/tau);
end

C = 1; % overall proportionality constant
rhoT = (-1i*C*Fsum)./(1+Fsum/tau);

plot(om,imag(rhoT));

xlabel('frequency');
ylabel('intensity');
