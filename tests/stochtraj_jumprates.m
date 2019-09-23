function [err,data] = test(opt,olddata)
% Check that using stochtraj_jump generates proper state residence times

Par.nTraj = 500;

rng_(1);

kp = rand()*1e9; % rate constant for forward process A -> B
km = rand()*1e9; % rate constant for reverse process B -> A
Sys.TransRates = [-kp, +km; +kp, -km].';

Sys.Orientations = [0,  0, 0;
                    0, pi, 0];

tauL = -1/mean(diag(Sys.TransRates));  % mean lifetime

D = eig(Sys.TransRates);
tau = -1/D(abs(D)/max(abs(D))>1e-11);  % relaxation time
            
Par.dt = tau/5;
Par.nSteps = ceil(200*tauL/Par.dt);

nSteps = Par.nSteps;

Opt.statesOnly = true;
[t,stateTraj] = stochtraj_jump(Sys,Par,Opt);

AutoCorrFFT = runprivate('autocorrfft', stateTraj.', 2, 1, 1).';

N = round(nSteps/2);

AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(N:end));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-t(:)/tau);

residuals = AutoCorrFFT - analytic;
rmsd = sqrt(mean(residuals(1:N).^2));

err = rmsd > 1e-2 || isnan(rmsd);

if opt.Display
  x = t(1:N)/1e-6;
  plot(x, AutoCorrFFT(1:N), x, analytic(1:N));
  xlabel('time (\mus)');
  legend('autocorr','analytical');
  axis tight
end

data = [];

end
