function ok = test(opt)

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

c20 = 3;
Sys.Potential = [2, 0, 0, c20];
[~, RTraj, qTraj] = stochtraj_diffusion(Sys,Par,Opt);

VecTraj = squeeze(RTraj(:, 3, :, :));

bins = linspace(0, pi, nBins)';
ThetaHist = zeros(nBins, nTraj);

N = round(nSteps/2);

for iTraj = 1:nTraj
  [alpha,beta,gamma] = quat2euler(qTraj(:,N:end,iTraj),'active');
  ThetaHist(:,iTraj) = hist(squeeze(beta), bins);
end

ThetaHist = mean(ThetaHist, 2);
ThetaHist = ThetaHist/sum(ThetaHist);

BoltzDist = exp(c20*(1.5*cos(bins).^2 - 0.5));
BoltzInt = sum(BoltzDist.*sin(bins));
BoltzDist = BoltzDist.*sin(bins)./BoltzInt;

ok = areequal(ThetaHist,BoltzDist,0.15,'rel');

if opt.Display
  plot(bins, ThetaHist, 'ro', bins, BoltzDist, 'k')
end
