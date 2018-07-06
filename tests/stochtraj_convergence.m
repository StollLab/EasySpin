function [err,data] = test(opt,olddata)
% Check that using stochtraj_diffusion with with convergence-checking works
% properly

Sys.tcorr = 10e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(50*Sys.tcorr/Par.dt);
Par.nTraj = 500;
Opt.chkcon = 1;
% Opt.Verbosity = 1;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

nBins = 50;

c20 = 3;
Sys.Potential.lambda = [c20, c20];
Sys.Potential.LMK = [2, 0, 0];
[~, RTraj] = stochtraj_diffusion(Sys,Par,Opt);

VecTraj = squeeze(RTraj(:, 3, :, :));

bins = linspace(0, pi, nBins)';
ThetaHist = zeros(nBins, nTraj);

N = round(nSteps/2);

for iTraj = 1:nTraj
  ThetaHist(:, iTraj) = hist(squeeze(acos(VecTraj(3,iTraj,N:end))), bins);
end

ThetaHist = mean(ThetaHist, 2);
ThetaHist = ThetaHist/sum(ThetaHist);

BoltzDist = exp(c20*(1.5*cos(bins).^2 - 0.5));
BoltzInt = sum(BoltzDist.*sin(bins));
BoltzDist = BoltzDist.*sin(bins)./BoltzInt;

residuals = ThetaHist - BoltzDist;
rmsd = sqrt(mean(residuals.^2));


if rmsd > 1e-2
  err = 1;
  plot(bins, ThetaHist, bins, BoltzDist)
else  
  err = 0;
end

data = [];

end
