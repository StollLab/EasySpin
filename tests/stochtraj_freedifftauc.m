function ok = test(opt)

% Check that using stochtraj_diffusion with free diffusion generates a proper
% rotational correlation time.

rng(1);

Sys.tcorr = 10e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;

[t, RTraj] = stochtraj_diffusion(Sys,Par);

VecTraj = squeeze(RTraj(:,3,:,:));

AutoCorrFFT = runprivate('autocorrfft', VecTraj.^2, 2, 1, 1);
AutoCorrFFT = mean(AutoCorrFFT, 3);


analytic = exp(-(1/Sys.tcorr)*t);
AutoCorrFFT = rescaledata(AutoCorrFFT,analytic,'lsq0');
AutoCorrFFT = AutoCorrFFT(:);

N = round(Par.nSteps/2);
residuals = AutoCorrFFT(1:N) - analytic(1:N);
rmsd = sqrt(mean(residuals.^2));

ok = rmsd<1e-2;

if opt.Display
  plot(t(1:N)/1e-6, AutoCorrFFT(1:N), t(1:N)/1e-6, analytic(1:N))
  xlabel('time (\mus)');
  legend('autocorrfft','analytical');
  axis tight
end
