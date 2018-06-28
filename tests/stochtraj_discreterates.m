function [err,data] = test(opt,olddata)
% Check that using stochtraj with Discrete model generates proper
% state residence times

Par.nTraj = 500;
Sys.TransRates = 1e9*[-0.5,  0.5;
                  0.5, -0.5];
Sys.States = [0,  0;
              0, pi;
              0,  0];

tau = -1/(2*Sys.TransRates(1,1));  % mean residence time
            
Par.dt = tau/5;
Par.nSteps = ceil(200*tau/Par.dt);

Opt.Model = 'Discrete';

TransRates = Sys.TransRates;
States = euler2quat(Sys.States);
nTraj = Par.nTraj;
nSteps = Par.nSteps;

[t, qTraj] = stochtraj(Sys,Par,Opt);
RTraj = quat2rotmat(qTraj);

VecTraj = squeeze(RTraj(:,3,:,:));

AutoCorrFFT = runprivate('autocorrfft', VecTraj.^2, 3, 1, 1);

AutoCorrFFT = squeeze(mean(AutoCorrFFT, 2));

N = round(nSteps/2);

AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(N:end));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-t/tau);

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
