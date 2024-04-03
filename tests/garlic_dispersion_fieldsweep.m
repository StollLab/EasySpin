function ok = test(opt)

% Calculate absorption and dispersion lineshapes using the solution
% of the Bloch equations
gam = gammae;  % gyromagnetic ratio, rad/s/T
B0 = 0.340 + linspace(-1,1,1e4)*0.002;  % static field, T
om = 2*pi*9.525*1e9;  % angular mw frequency, rad/s

T1 = 1e-6;  % s
T2 = 1e-7;  % s
B1 = 1e-6;  % T
omL = -gam*B0;

% Complex steady-state solution of Bloch equations
chi = gam*B0*T2.*((om-omL)*T2+1i)./(1+(om-omL).^2*T2^2+gam^2*B1^2*T1*T2);
chi1 = real(chi);  % dispersion
chi2 = -imag(chi);  % absorption
chi1 = chi1/max(abs(chi1));
chi2 = chi2/max(abs(chi2));

% Calculate lineshape using pepper
Exp.mwFreq = om/2/pi/1e9;
Exp.Range = [min(B0) max(B0)]*1e3;
Exp.nPoints = numel(B0);
Exp.Harmonic = 0;

Sys.g = gfree;
fwhm = 1/pi/T2; % Hz
fwhm = fwhm*planck/bmagn/Sys.g/1e-3; % Hz -> mT
Sys.lw = [0 fwhm];  % mT
[~,chi2_sim] = garlic(Sys,Exp);
chi2_sim = chi2_sim/max(abs(chi2_sim));

Exp.mwPhase = pi/2;
[~,chi1_sim] = garlic(Sys,Exp);
chi1_sim = chi1_sim/max(abs(chi1_sim));

if opt.Display
  tiledlayout(2,1)
  nexttile
  plot(B0,chi1,B0,chi1_sim,'--',B0,(chi1_sim-chi1)*10);
  axis tight
  xlabel('field (T)')
  legend('explicit','pepper','difference*10')
  title('dispersion')

  nexttile
  plot(B0,chi2,B0,chi2_sim,'--',B0,(chi2_sim-chi2)*10);
  axis tight
  xlabel('field (T)')
  legend('explicit','pepper','difference*10')
  title('absorption')
end

threshold = 4e-2;
ok(1) = areequal(chi1,chi1_sim,threshold,'rel');
ok(2) = areequal(chi2,chi2_sim,threshold,'rel');
