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
  AutoCorrFFT(iTraj,:) = autocorrfft(squeeze(VecTraj(:,iTraj,:).^2));
%   AutoCorrFFT(:, iTraj) = autocorrfft(q(:, :, iTraj).^2);
end

AutoCorrFFT = sum(AutoCorrFFT, 1)'/nTraj;
AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(end/2:end));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-(1/tcorr)*t);

% ChiSquare = sum(((AutoCorrFFT - analytic).^2)./analytic)
rmsd = sqrt(sum((AutoCorrFFT - analytic).^2)/nSteps);

if rmsd > 1e-2
  err = 1;
  plot(t, AutoCorrFFT, t, analytic)
else  
  err = 0;
end

data = [];

end
