function [err,data] = test(opt,olddata)
% Check that using stochtraj with free diffusion generates a proper 
% distribution of orientations

Par.tcorr = 10*rand()*1e-9;
Par.dt = Par.tcorr/10;
Par.nSteps = ceil(100*Par.tcorr/Par.dt);
Par.nTraj = 100;
Par.theta = pi*(2*rand()-1);
Par.phi = 2*pi*(2*rand()-1);
Par.chi = 2*pi*(2*rand()-1);

nTraj = Par.nTraj;
nSteps = Par.nSteps;
c20 = 0.0;

nBins = 70;

[t, R] = stochtraj(Par);


VecTraj = squeeze(R(:, 3, :, :));

bins = linspace(0, pi, nBins)';
ThetaHist = zeros(nBins, nTraj);

for iTraj = 1:nTraj
  ThetaHist(:, iTraj) = hist(acos(VecTraj(3, :, iTraj)), bins);
end

% ThetaHistAvg = sum(ThetaHist, 2)/nTraj;
% HistDiff = bsxfun(@minus,ThetaHist,ThetaHistAvg);
% rmsd = sqrt(sum(HistDiff.^2,1)/70);
% mean(rmsd)

% plot(hist(rmsd, 70))
% plot(rmsd)
ThetaHist = ThetaHist/sum(ThetaHist);

BoltzDist = exp(c20*(1.5*cos(bins).^2 - 0.5));
BoltzInt = sum(BoltzDist.*sin(bins));
BoltzDist = BoltzDist.*sin(bins)./BoltzInt;

%ChiSquare = sum(((ThetaHist - BoltzDist).^2)./ThetaHist)
rmsd = sqrt(sum((ThetaHist - BoltzDist).^2)/nBins);

if rmsd > 1e-2
  err = 1;
  plot(bins, ThetaHist, bins, BoltzDist)
else  
  err = 0;
end

data = [];

end
