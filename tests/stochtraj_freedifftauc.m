function [err,data] = test(opt,olddata)
% Check that using stochtraj with free diffusion generates a proper
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

[t, R, q] = stochtraj(Sys,Par);

VecTraj = squeeze(R(:,3,:,:));

AutoCorrFFT = zeros(nTraj,nSteps);

for iTraj = 1:nTraj
  AutoCorrFFT(iTraj,:) = autocorrfft(squeeze(VecTraj(:,iTraj,:).^2), 1);
%   AutoCorrFFT(:, iTraj) = autocorrfft(q(:, :, iTraj).^2);
end

N = round(nSteps/2);

AutoCorrFFT = mean(AutoCorrFFT, 1).';
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
