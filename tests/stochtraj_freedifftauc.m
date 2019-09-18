function [err,data] = test(opt,olddata)
% Check that using stochtraj_diffusion with free diffusion generates a proper
% rotational correlation time

Sys.tcorr = 10*rand()*1e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;
Par.OriStart = pi*[1 2 2].*(2*rand(1,3)-1);

tcorr = Sys.tcorr;
nSteps = Par.nSteps;

[t, RTraj] = stochtraj_diffusion(Sys,Par);

VecTraj = squeeze(RTraj(:,3,:,:));

AutoCorrFFT = runprivate('autocorrfft', VecTraj.^2, 2, 1, 1);

AutoCorrFFT = mean(AutoCorrFFT, 3);

N = round(nSteps/2);

AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(N:end));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-(1/tcorr)*t);

residuals = AutoCorrFFT - analytic;
rmsd = sqrt(mean(residuals(1:N).^2));

err = rmsd>1e-2;

if opt.Display
  plot(t(1:N)/1e-6, AutoCorrFFT(1:N), t(1:N)/1e-6, analytic(1:N))
  xlabel('time (\mus)');
  legend('autocorrfft','analytical');
end

data = [];

end
