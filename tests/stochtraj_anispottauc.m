function [err,data] = test(opt,olddata)
% Check that using stochtraj_diffusion with anisotropic diffusion generates a
% proper rotational correlation timez

Sys.tcorr = 10*rand()*1e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;
Par.Omega = [  pi*(2*rand()-1); 
             2*pi*(2*rand()-1);
             2*pi*(2*rand()-1) ];

tcorr = Sys.tcorr;
nTraj = Par.nTraj;
nSteps = Par.nSteps;

c20 = 10;
Sys.Potential = [2, 0, 0, c20];
[t, RTraj] = stochtraj_diffusion(Sys,Par);

VecTraj = squeeze(RTraj(:,3,:,:));

AutoCorrFFT = runprivate('autocorrfft', VecTraj.^2, 3, 1, 1);

AutoCorrFFT = squeeze(mean(AutoCorrFFT, 2));

N = round(nSteps/2);
M = round(N/2);

AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(M:N));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-(1/tcorr)*t);

% ChiSquare = sum(((AutoCorrFFT - analytic).^2)./AutoCorrFFT)
residuals = AutoCorrFFT - analytic;
rmsd = sqrt(mean(residuals(1:N).^2));

[k, c, yFit] = exponfit(t(1:N), AutoCorrFFT(1:N));
tauR = 1/k;

residuals = AutoCorrFFT(1:N) - yFit;
rmsd = sqrt(mean(residuals.^2));

if rmsd > 1e-2 || isnan(rmsd) || tcorr-tauR < 0
  err = 1;
  plot(t(1:N)/1e-6, AutoCorrFFT(1:N), t(1:N)/1e-6, analytic(1:N))
  xlabel('t (\mu s)')
  subplot(2,1,2)
  plot(t(1:N)/1e-6, cumtrapz(t(1:N), AutoCorrFFT(1:N)))
  xlabel('t (\mu s)')
else
  err = 0;
end

data = [];

end
