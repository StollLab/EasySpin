function [err,data] = test(opt,olddata)
% Check that using stochtraj_diffusion with with convergence-checking works
% properly

Sys.tcorr = 10e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 500;
Opt.chkcon = 1;
% Opt.Verbosity = 1;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

nBins = 50;

N = 100;
alphaGrid = linspace(0, 2*pi, N);
betaGrid = linspace(0, pi, N);
gammaGrid = linspace(0, 2*pi, N);

[AlphaGrid,BetaGrid,GammaGrid] = ndgrid(alphaGrid,betaGrid,gammaGrid);

c20 = 3;
LMK = [2,0,0];
Yfun = real(wignerd(LMK,AlphaGrid,BetaGrid,GammaGrid));
U = -c20*Yfun;

Sys.Potential = U;

% c20 = 3;
% Sys.Potential = [2, 0, 0, c20];
[~, RTraj, qTraj] = stochtraj_diffusion(Sys,Par,Opt);

VecTraj = squeeze(RTraj(:, 3, :, :));

bins = linspace(0, pi, nBins)';
ThetaHist = zeros(nBins, nTraj);

N = round(nSteps/2);

for iTraj = 1:nTraj
  [alpha,beta,gamma] = quat2euler(qTraj(:,iTraj,N:end));
  ThetaHist(:,iTraj) = hist(squeeze(beta), bins);
%   ThetaHist(:, iTraj) = hist(squeeze(acos(VecTraj(3,iTraj,N:end))), bins);
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
