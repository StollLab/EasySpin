function ok = test(opt)

% Calculate absorption and dispersion lineshapes using the solution
% of the Bloch equations
gam = gammae;  % gyromagnetic ratio, rad/s/T
B0 = 0.340;  % static field, T
om = 2*pi*linspace(9.5,9.56,10001)*1e9;  % angular mw frequency, rad/s

T1 = 1e-6;  % s
T2 = 1e-7;  % s
B1 = 1e-6;  % T
omL = -gam*B0;

% Complex steady-state solution of Bloch equations
chi = gam*B0*T2*((om-omL)*T2+1i)./(1+(om-omL).^2*T2^2+gam^2*B1^2*T1*T2);
chi1 = real(chi);  % dispersion
chi2 = -imag(chi);  % absorption
chi1 = chi1/max(abs(chi1));
chi2 = chi2/max(abs(chi2));

% Calculate lineshape using pepper
nu_GHz = om/2/pi/1e9;
Exp.mwRange = [min(nu_GHz) max(nu_GHz)];
Exp.nPoints = numel(nu_GHz);
Exp.Field = B0*1e3;

Sys.g = gfree;
fwhm = 1/pi/T2;
Sys.lw = [0.01 fwhm/1e6];  % MHz; small non-zero Gaussian broadening needed
[~,chi2_sim] = pepper(Sys,Exp);
chi2_sim = chi2_sim/max(abs(chi2_sim));

Exp.mwPhase = pi/2;
[nu,chi1_sim] = pepper(Sys,Exp);
chi1_sim = chi1_sim/max(abs(chi1_sim));

if opt.Display
  tiledlayout(2,1)
  nexttile
  plot(nu_GHz,chi1,nu_GHz,chi1_sim,'--');
  axis tight
  xlabel('frequency (GHz)')
  legend('explicit','pepper')
  title('dispersion')

  nexttile
  plot(nu_GHz,chi2,nu_GHz,chi2_sim,'--');
  axis tight
  xlabel('frequency (GHz)')
  legend('explicit','pepper')
  title('absorption')
end

threshold = 1e-3;
ok(1) = areequal(chi1,chi1_sim,threshold,'rel');
ok(2) = areequal(chi2,chi2_sim,threshold,'rel');
