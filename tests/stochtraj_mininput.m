function ok = test(opt)

% Check that using stochtraj_diffusion with free diffusion generates a proper 
% distribution of orientations

Sys.tcorr = 10*rand()*1e-9;
Par.nTraj = 100;

nTraj = Par.nTraj;

nBins = 70;

[t, RTraj] = stochtraj_diffusion(Sys,Par);
nSteps = size(RTraj,3);

VecTraj = squeeze(RTraj(:, 3, :, :));

bins = linspace(0, pi, nBins)';
ThetaHist = zeros(nBins, nTraj);

for iTraj = 1:nTraj
  ThetaHist(:, iTraj) = hist(acos(VecTraj(3, :, iTraj)), bins);
end

ThetaHist = ThetaHist/sum(ThetaHist);

BoltzDist = ones(nBins,1);
BoltzInt = sum(BoltzDist.*sin(bins));
BoltzDist = BoltzDist.*sin(bins)./BoltzInt;

%ChiSquare = sum(((ThetaHist - BoltzDist).^2)./ThetaHist)
residuals = ThetaHist - BoltzDist;
rmsd = sqrt(mean(residuals.^2));

ok = rmsd < 1e-2;

if opt.Display
  plot(bins, ThetaHist, bins, BoltzDist)
end
