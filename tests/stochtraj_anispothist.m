function [err,data] = test(opt,olddata)

Par.lambda = 2*rand();
Par.tcorr = 10*rand()*1e-9;
Par.dt = 0.5e-9;
Par.nSteps = 4000;
Par.nTraj = 800;

nTraj = Par.nTraj;
nSteps = Par.nSteps;
c20 = Par.lambda;

nBins = 50;

[t, R] = stochtraj(Par);


VecTraj = squeeze(R(:, 3, :, :));

bins = linspace(0, pi, nBins)';
ThetaHist = zeros(nBins, nTraj);

for iTraj = 1:nTraj
  ThetaHist(:, iTraj) = hist(acos(VecTraj(3, :, iTraj)), bins);
end

ThetaHist = sum(ThetaHist, 2);
ThetaHist = ThetaHist/sum(ThetaHist);

BoltzDist = exp(c20*(1.5*cos(bins).^2 - 0.5));
BoltzInt = sum(BoltzDist.*sin(bins));
BoltzDist = BoltzDist.*sin(bins)./BoltzInt;

ChiSquare = sum(((ThetaHist - BoltzDist).^2)./ThetaHist);

% This seems like a loose condition and should be investigated further
if ChiSquare > 5e-3
  err = 1;
else  
  err = 0;
end

data = [];

end
