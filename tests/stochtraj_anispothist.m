function [err,data] = test(opt,olddata)
% Check that using stochtraj with anisotropic diffusion generates a
% proper distribution of orientations

Sys.tcorr = 10e-9;%10*rand()*1e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(300*Sys.tcorr/Par.dt);
Par.nTraj = 400;
Par.beta = 0;%pi*(2*rand()-1);
Par.alpha = 0;%2*pi*(2*rand()-1);
Par.gamma = 0;%2*pi*(2*rand()-1);
% [Par.alpha, Par.beta, Par.gamma]/pi*180

nTraj = Par.nTraj;
nSteps = Par.nSteps;

nBins = 50;

% c20 = 2*(2*rand()-1);
% Sys.lambda = struct('c20', c20);
% [t, R_c20] = stochtraj(Sys,Par);
% Sys.lambda = rmfield(Sys.lambda,'c20');

% c2p1 = 2*(2*rand()-1);
% Sys.lambda = struct('c2p1', c2p1);
% [t, R_c2p1] = stochtraj(Sys,Par);
% Sys.lambda = rmfield(Sys.lambda,'c2p1');

c2p2 = -1;%2*(2*rand()-1)
Sys.lambda = struct('c2p2', c2p2, 'c2m2', c2p2);
[t, R_c2p2] = stochtraj(Sys,Par);
Sys.lambda = rmfield(Sys.lambda,'c2p2');
Sys.lambda = rmfield(Sys.lambda,'c2m2');


% VecTraj_c20 = squeeze(R_c20(:,3,:,:));
% VecTraj_c2p1 = squeeze(R_c2p1(:,3,:,:));
VecTraj_c2p2 = squeeze(R_c2p2(:,3,:,:));

bins = linspace(0, pi, nBins)';
% ThetaHist_c20 = zeros(nBins, nTraj);
% ThetaHist_c2p1 = zeros(nBins, nTraj);
ThetaHist_c2p2 = zeros(nBins, nTraj);

for iTraj = 1:nTraj
%   ThetaHist_c20(:,iTraj) = hist(squeeze(acos(VecTraj_c20(3,iTraj,:))), bins);
%   ThetaHist_c2p1(:,iTraj) = hist(squeeze(acos(VecTraj_c2p1(3,iTraj,:))), bins);
  ThetaHist_c2p2(:,iTraj) = hist(squeeze(acos(VecTraj_c2p2(3,iTraj,:))), bins);
end

% ThetaHist_c20 = sum(ThetaHist_c20, 2);
% ThetaHist_c2p1 = sum(ThetaHist_c2p1, 2);
ThetaHist_c2p2 = sum(ThetaHist_c2p2, 2);

% ThetaHist_c20 = ThetaHist_c20/sum(ThetaHist_c20);
% ThetaHist_c2p1 = ThetaHist_c2p1/sum(ThetaHist_c2p1);
ThetaHist_c2p2 = ThetaHist_c2p2/sum(ThetaHist_c2p2);

% BoltzDist_c20 = exp(c20*1/2*(3*cos(bins).^2 - 1));
% BoltzInt_c20 = sum(BoltzDist_c20.*sin(bins));
% BoltzDist_c20 = BoltzDist_c20.*sin(bins)./BoltzInt_c20;

% BoltzDist_c2p1 = exp(c2p1*sqrt(3/8)*sin(2*bins));
% BoltzInt_c2p1 = sum(BoltzDist_c2p1.*sin(bins));
% BoltzDist_c2p1 = BoltzDist_c2p1.*sin(bins)./BoltzInt_c2p1;

BoltzDist_c2p2 = exp(c2p2*2*sqrt(3/8)*sin(bins).^2);
BoltzInt_c2p2 = sum(BoltzDist_c2p2.*sin(bins));
BoltzDist_c2p2 = BoltzDist_c2p2.*sin(bins)./BoltzInt_c2p2;

% ChiSquare = sum(((ThetaHist - BoltzDist).^2)./ThetaHist);
% rmsd_c20 = sqrt(sum((ThetaHist_c20 - BoltzDist_c20).^2)/nBins);
% rmsd_c2p1 = sqrt(sum((ThetaHist_c2p1 - BoltzDist_c2p1).^2)/nBins);
rmsd_c2p2 = sqrt(sum((ThetaHist_c2p2 - BoltzDist_c2p2).^2)/nBins);

% This seems like a loose condition and should be investigated further
if rmsd_c2p2 > 1e-2% || rmsd_c2p1 > 1e-2 || rmsd_c2p2 > 1e-2
  err = 1;
  hold on
%   plot(bins, ThetaHist_c20, bins, BoltzDist_c20)
%   plot(bins, ThetaHist_c2p1, bins, BoltzDist_c2p1)
%   plot(bins, ThetaHist_c2p2, bins, BoltzDist_c2p2)
  hold off
else  
  err = 0;
  plot(bins, ThetaHist_c2p2, bins, BoltzDist_c2p2)
end

data = [];

end
