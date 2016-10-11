function [err,data] = test(opt,olddata)
% Check that using stochtraj with free diffusion generates a proper
% rotational correlation time

Par.tcorr = 10*rand()*1e-9;
Par.dt = Par.tcorr/10;
Par.nSteps = ceil(100*Par.tcorr/Par.dt);
Par.nTraj = 800;
Par.theta = pi*(2*rand()-1);
Par.phi = 2*pi*(2*rand()-1);
Par.chi = 2*pi*(2*rand()-1);

tcorr = Par.tcorr;
nTraj = Par.nTraj;
nSteps = Par.nSteps;

[t, R, q] = stochtraj(Par);

VecTraj = squeeze(R(:, 3, :, :));

AutoCorrFFT = zeros(nSteps, nTraj);

for iTraj = 1:nTraj
  AutoCorrFFT(:, iTraj) = autocorrfft(VecTraj(:, :, iTraj).^2);
%   AutoCorrFFT(:, iTraj) = autocorrfft(q(:, :, iTraj).^2);
end

AutoCorrFFT = sum(AutoCorrFFT, 2)'/nTraj;

analytic = exp(-(1/tcorr)*t);

% ChiSquare = sum(((AutoCorrFFT - analytic).^2)./analytic)
rmsd = sqrt(sum((AutoCorrFFT - analytic).^2)/nSteps);

if rmsd > 5e-2
  err = 1;
  plot(t, AutoCorrFFT, t, analytic)
else  
  err = 0;
end

data = [];

end
