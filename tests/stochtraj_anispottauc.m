function [err,data] = test(opt,olddata)
% Check that using stochtraj with anisotropic diffusion generates a
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

c20 = 2*rand();
Sys.Coefs = [c20, c20];
Sys.LMK = [2, 0, 0];
[t, R] = stochtraj(Sys,Par);

VecTraj = squeeze(R(:,3,:,:));

AutoCorrFFT = zeros(nTraj, nSteps);

for iTraj = 1:nTraj
  AutoCorrFFT(iTraj,:) = autocorrfft(squeeze(VecTraj(:,iTraj,:).^2), 1);
end

N = round(nSteps/2);
M = round(N/2);

AutoCorrFFT = mean(AutoCorrFFT, 1).';
AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(M:N));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-(1/tcorr)*t);

% ChiSquare = sum(((AutoCorrFFT - analytic).^2)./AutoCorrFFT)
residuals = AutoCorrFFT - analytic;
rmsd = sqrt(mean(residuals(1:N).^2));

if rmsd > 1e-2 || isnan(rmsd)
  err = 1;
%   subplot(2,1,1)
%   plot(t(1:N)/1e-6, AutoCorrFFT(1:N), t(1:N)/1e-6, analytic(1:N))
%   xlabel('t (\mu s)')
%   subplot(2,1,2)
%   plot(t(1:N)/1e-6, cumtrapz(t(1:N), AutoCorrFFT(1:N)))
%   xlabel('t (\mu s)')
else
  err = 0;  

end

data = [];

end
