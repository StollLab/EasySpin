function [err,data] = test(opt,olddata)
% Check that using stochtraj with anisotropic diffusion generates a
% proper distribution of orientations

Sys.tcorr = 10e-9;
Par.dt = Sys.tcorr/5;
Par.nSteps = ceil(50*Sys.tcorr/Par.dt);
Par.nTraj = 200;
Par.beta = pi*(2*rand()-1);
Par.alpha = 2*pi*(2*rand()-1);
Par.gamma = 2*pi*(2*rand()-1);
Par.chkcon = 1;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

nBins = 50;

c20 = 3;
Sys.Coefs = [c20, c20];
Sys.LMK = [2, 0, 0];
[t, R] = stochtraj(Sys,Par);


VecTraj = squeeze(R(:, 3, :, :));

bins = linspace(0, pi, nBins)';
ThetaHist = zeros(nBins, nTraj);

for iTraj = 1:nTraj
  ThetaHist(:, iTraj) = hist(squeeze(acos(VecTraj(3, iTraj, :))), bins);
end

ThetaHist = sum(ThetaHist, 2);
ThetaHist = ThetaHist/sum(ThetaHist);

BoltzDist = exp(c20*(1.5*cos(bins).^2 - 0.5));
BoltzInt = sum(BoltzDist.*sin(bins));
BoltzDist = BoltzDist.*sin(bins)./BoltzInt;

% ChiSquare = sum(((ThetaHist - BoltzDist).^2)./ThetaHist);
rmsd = sqrt(sum((ThetaHist - BoltzDist).^2)/nBins);

% This seems like a loose condition and should be investigated further
if rmsd > 1e-2
  err = 1;
  plot(bins, ThetaHist, bins, BoltzDist)
else  
  err = 0;
end

data = [];

end
