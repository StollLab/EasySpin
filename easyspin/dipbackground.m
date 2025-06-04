% dipbackground  Intermolecular dipolar background decay
%
%  Vinter = dipbackground(t,conc,lambda)
%  dipbackground(...)
%
%  Calculates the intermolecular dipolar background decay over the time axis,
%  for spin concentration conc (in µM) and modulation depth lambda (between 0 and 1).
%
%  Input:
%    t       time axis (µs)
%    conc    spin concentration (µM)
%    lambda  modulation depth
%
%  Output:
%    Vinter  intermolecular background decay function
%
%  Example:
%    t = linspace(-0.2,5);  % µs
%    Vinter = dipbackground(t,100,0.3);
%    plot(t,Vinter);

function Vinter = dipbackground(t,conc_uM,lambda)

if nargin==0
  help(mfilename);
  return
end

g1 = gfree;
g2 = gfree;

D = mu0/(4*pi)*g1*g2*bmagn^2/hbar;  % rad/s
conc_mM = conc_uM/1e3;  % µM -> mM=mol/m^3
conc = conc_mM*avogadro;  % mM=mol/m^3 -> m^-3
t = t*1e-6;  % µs -> s
k = 8*pi^2/(9*sqrt(3))*D*conc*lambda;
Vinter = exp(-k*abs(t));

% Plotting
if nargout==0
  plot(t/1e-6,Vinter)
  xlabel('time (µs)')
  ylabel('{\it{V}}_{inter}')
  grid on
end

end
