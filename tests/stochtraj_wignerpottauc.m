function [err,data] = test(opt,olddata)
% Check that using stochtraj_diffusion with anisotropic diffusion generates a
% proper rotational correlation timez

rng(1)

Sys.tcorr = 10*1e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;
Par.OriStart = [  pi*rand(); 
                2*pi*rand();
                2*pi*rand()];

tcorr = Sys.tcorr;
nTraj = Par.nTraj;
nSteps = Par.nSteps;

c20 = 2;
Sys.Potential = [2, 0, 0, c20];
[t, RTraj] = stochtraj_diffusion(Sys,Par);

VecTraj = squeeze(RTraj(:,3,:,:));

AutoCorrFFT = runprivate('autocorrfft', VecTraj.^2, 3, 1, 1);

AutoCorrFFT = squeeze(mean(AutoCorrFFT, 2));

N = round(nSteps/2);
M = round(N/2);

AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(M:N));
y = AutoCorrFFT/max(AutoCorrFFT);

data.y = y;

if ~isempty(olddata)
  err = any(abs(olddata.y-y)>1e-10);
else
  err = [];
end

end
