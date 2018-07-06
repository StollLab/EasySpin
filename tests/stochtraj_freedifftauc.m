function [err,data] = test(opt,olddata)
% Check that using stochtraj_diffusion with free diffusion generates a proper
% rotational correlation time

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

[t, RTraj] = stochtraj_diffusion(Sys,Par);

VecTraj = squeeze(RTraj(:,3,:,:));

AutoCorrFFT = runprivate('autocorrfft', VecTraj.^2, 3, 1, 1);

AutoCorrFFT = squeeze(mean(AutoCorrFFT, 2));

N = round(nSteps/2);

AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(N:end));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-(1/tcorr)*t);

residuals = AutoCorrFFT - analytic;
rmsd = sqrt(mean(residuals(1:N).^2));

if rmsd > 1e-2
  err = 1;
  plot(t(1:N), AutoCorrFFT(1:N), t(1:N), analytic(1:N))
else  
  err = 0;
end

data = [];

end
