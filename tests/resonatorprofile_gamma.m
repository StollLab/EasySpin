function ok = test(opt)

% Test voltage reflection against expressions from
%   Vladimir Krymov and Gary Gerfen
%   Analysis of the tuning and operation of reflection resonator EPR spectrometers
%   J. Magn. Reson. 162 (2003) 4660478
%   https://doi.org/10.1016%2FS1090-7807(03)00109-5

% Settings
nu = linspace(9,10,201);  % GHz
nu0 = 9.6;  % GHz
Qu = 200;
beta = 3.2;

% Expressions from the paper 
xi = -Qu*(nu/nu0-nu0./nu);  % Eq.(6)
Gamma = (1-beta-1i*xi)./(1+beta-1i*xi);  % Eq.(7)

Gamma_ = resonatorprofile(nu,nu0,Qu,beta,'voltagereflection');

ok = areequal(Gamma,Gamma_,1e-12,'abs');

if opt.Display
  plot(nu,real(Gamma),nu,real(Gamma_),'.',nu,imag(Gamma),nu,imag(Gamma_),'.');
  xlabel('frequency (GHz)');
  ylabel('Re(\Gamma), Im(\Gamma)');
  legend('Re(\Gamma) ref','Re(\Gamma)','Im(\Gamma) ref','Im(\Gamma)','AutoUpdate','off');
  axis tight
  grid on
  xline(nu0);
end
