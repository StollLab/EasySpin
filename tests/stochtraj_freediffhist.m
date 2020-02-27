function ok = test(opt)

% Check that using stochtraj_diffusion with free diffusion generates a proper 
% distribution of orientations

rng(1);
Sys.tcorr = 10*rand()*1e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;
Par.OriStart = [  pi*(2*rand()-1); 
                2*pi*(2*rand()-1);
                2*pi*(2*rand()-1) ];

nTraj = Par.nTraj;
nSteps = Par.nSteps;

nBins = 50;

[t, RTraj] = stochtraj_diffusion(Sys,Par);

VecTraj = squeeze(RTraj(:, 3, :, :));

bins = linspace(0, pi, nBins)';
ThetaHist = zeros(nBins, nTraj);

for iTraj = 1:nTraj
  ThetaHist(:, iTraj) = hist(squeeze(acos(VecTraj(3, iTraj, :))), bins);
end

ThetaHist = mean(ThetaHist, 2);
ThetaHist = ThetaHist/sum(ThetaHist);

BoltzDist = ones(nBins,1);
BoltzInt = sum(BoltzDist.*sin(bins));
BoltzDist = BoltzDist.*sin(bins)./BoltzInt;

residuals = ThetaHist - BoltzDist;
rmsd = sqrt(mean(residuals.^2));

ok = rmsd<1e-2;

if opt.Display
  plot(bins, ThetaHist, bins, BoltzDist)
end
