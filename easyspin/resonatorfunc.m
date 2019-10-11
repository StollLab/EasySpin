% resonatorfunc  Resonator reflection coefficient
%
%  resonatorfunc(nu,nu0,Qu,beta)
%  Gamma = resonatorfunc(nu,nu0,Qu,beta)
%
%  Calculates the reflection coefficient for a reflection resonator, i.e. the
%  fraction of voltage from an incoming wave that is reflected by the resonator.
%
%  Inputs:
%    nu      frequency range array, GHz
%    nu0     resonator frequency, GHz
%    Qu      unloaded Q-factor of the resonator
%    beta    coupling coefficient
%              beta<1: undercoupled
%              beta=1: critically coupled (matched)
%              beta>1: overcoupled
%
%    If nu is [], then an appropriate frequency range is chosen automatically.
%
%  Outputs:
%    Gamma   Voltage reflection coefficient (ratio of E-field amplitudes
%            of reflected and incoming wave). Gamma is the _voltage_ reflection
%            coefficient. The power reflection coefficient is abs(Gamma)^2.
%
%    If no output is requested, the results are plotted.
%
%  Example:
%    nu = linspace(9.2,9.8,1001);
%    nu0 = 9.5;
%    Qu = 1000;
%    beta = 1;
%    resonatorfunc(nu,nu0,Qu,beta);

function varargout = resonatorfunc(nu,nu0,Qu,beta)

if nargin==0
  help(mfilename);
  return
end

if isempty(nu)
  dnu = nu0/Qu*10;
  nu = linspace(max(nu0-dnu,0),nu0+dnu,1001);
end

xi = -Qu*(nu/nu0-nu0./nu);
Gamma = (1-beta-1i*xi)./(1+beta-1i*xi);  % voltage reflection coefficient

switch nargout
  case 0
    plot(nu,abs(Gamma));            % reflected voltage
    xlabel('frequency (GHz)');
    ylabel('|\Gamma|, norm of voltage reflection coefficient');
    title(sprintf('Qu = %d, beta = %g',Qu,beta));
    ylim([0 1]);
    line([1 1]*nu0,ylim,'Color','r');
  case 1
    varargout = {Gamma};
end
