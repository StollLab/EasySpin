function [err,data] = test(opt,olddata)
% Check that using stochtraj with anisotropic diffusion generates a
% proper rotational correlation time

Sys.lambda = 2*(2*rand()-1);
Sys.tcorr = 10*rand()*1e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(100*Sys.tcorr/Par.dt);
Par.nTraj = 800;
Par.theta = pi*rand();
Par.phi = 2*pi*rand();
Par.chi = 2*pi*rand();

tcorr = Sys.tcorr;
nTraj = Par.nTraj;
nSteps = Par.nSteps;

[t, R] = stochtraj(Sys,Par);

VecTraj = squeeze(R(:,3,:,:));

AutoCorrFFT = zeros(nSteps, nTraj);

for iTraj = 1:nTraj
  AutoCorrFFT(:, iTraj) = autocorrfft(VecTraj(:,:,iTraj).^2);
end

AutoCorrFFT = sum(AutoCorrFFT, 2).'/nTraj;

analytic = exp(-(1/tcorr)*t);

% ChiSquare = sum(((AutoCorrFFT - analytic).^2)./AutoCorrFFT)
rmsd = sqrt(sum((AutoCorrFFT - analytic).^2)/nSteps);

if rmsd > 5e-2
  err = 1;
  plot(t/1e-6, AutoCorrFFT, t/1e-6, analytic)
  xlabel('t (\mu s)')
else
  err = 0;
end

data = [];

end
