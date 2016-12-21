function [err,data] = test(opt,olddata)
% Check that using stochtraj with anisotropic diffusion generates a
% proper distribution of orientations

Sys.tcorr = 10*rand()*1e-9;;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;
Par.beta = pi*(2*rand()-1);
Par.alpha = 2*pi*(2*rand()-1);
Par.gamma = 2*pi*(2*rand()-1);

nTraj = Par.nTraj;
nSteps = Par.nSteps;

nBins = 50;

% c20 = 2;%2*(2*rand()-1);
% Sys.Coefs = [c20, c20];
% Sys.LMK = [2, 0, 0];
% [t, R_c20] = stochtraj(Sys,Par);

% c2p1 = 3;
% Sys.Coefs = [c2p1, c2p1];
% Sys.LMK = [2, 0, 1];
% [t, R_c2p1] = stochtraj(Sys,Par);

c2p2 = 2;
Sys.Coefs = [c2p2, c2p2];
Sys.LMK = [2, 0, 2];
[t, R_c2p2] = stochtraj(Sys,Par);


% plot(squeeze(atan2(VecTraj_c2p2(2,1,:),VecTraj_c2p2(1,1,:))))
% plot(squeeze(acos(VecTraj_c2p2(3,1,:)))/pi*180)

ThetaBins = linspace(0, pi, nBins)';
PhiBins = linspace(-pi, pi, nBins)';

% VecTraj_c20 = squeeze(R_c20(:,3,:,:));
% ThetaHist_c20 = zeros(nBins, nTraj);
% PhiHist_c20 = zeros(nBins, nTraj);

% VecTraj_c2p1 = squeeze(R_c2p1(:,3,:,:));
% ThetaHist_c2p1 = zeros(nBins, nTraj);
% PhiHist_c2p1 = zeros(nBins, nTraj);
% Hist2D_c2p1 = zeros(nBins-1, nBins-1, nTraj);

VecTraj_c2p2 = squeeze(R_c2p2(:,3,:,:));
ThetaHist_c2p2 = zeros(nBins, nTraj);
PhiHist_c2p2 = zeros(nBins, nTraj);
Hist2D_c2p2 = zeros(nBins-1, nBins-1, nTraj);

for iTraj = 1:nTraj
%   ThetaHist_c20(:,iTraj) = hist(squeeze(acos(VecTraj_c20(3,iTraj,:))),ThetaBins);
%   PhiHist_c20(:,iTraj) = hist(squeeze(atan2(VecTraj_c20(2,iTraj,:),...
%                                             VecTraj_c20(1,iTraj,:))),PhiBins);

%   theta = squeeze(acos(VecTraj_c2p1(3,iTraj,:)));
%   phi = squeeze(atan2(VecTraj_c2p1(2,iTraj,:), VecTraj_c2p1(1,iTraj,:)));
%   ThetaHist_c2p1(:,iTraj) = hist(theta, ThetaBins);
%   PhiHist_c2p1(:,iTraj) = hist(phi, PhiBins);
%   [Hist2D_c2p1(:,:,iTraj),~,~] = histcounts2(theta,phi,ThetaBins,PhiBins);
                                          
  theta = squeeze(acos(VecTraj_c2p2(3,iTraj,:)));
  phi = squeeze(atan2(VecTraj_c2p2(2,iTraj,:), VecTraj_c2p2(1,iTraj,:)));
  ThetaHist_c2p2(:,iTraj) = hist(squeeze(acos(VecTraj_c2p2(3,iTraj,:))), ThetaBins);
  PhiHist_c2p2(:,iTraj) = hist(squeeze(atan2(VecTraj_c2p2(2,iTraj,:),...
                                             VecTraj_c2p2(1,iTraj,:))), PhiBins);
  [Hist2D_c2p2(:,:,iTraj),~,~] = histcounts2(theta,phi,ThetaBins,PhiBins);

end

% ThetaHist_c20 = sum(ThetaHist_c20,2);
% PhiHist_c20 = sum(PhiHist_c20,2);
% ThetaHist_c20 = ThetaHist_c20/sum(ThetaHist_c20);
% PhiHist_c20 = PhiHist_c20/sum(PhiHist_c20);
% BoltzDist_c20 = exp(c20*1/2*(3*cos(ThetaBins).^2 - 1));
% BoltzInt_c20 = sum(BoltzDist_c20.*sin(ThetaBins));
% BoltzDist_c20 = BoltzDist_c20.*sin(ThetaBins)./BoltzInt_c20;
% rmsd_c20 = sqrt(sum((ThetaHist_c20 - BoltzDist_c20).^2)/nBins);

% Hist2D_c2p1 = sum(Hist2D_c2p1, 3);
% ThetaHist_c2p1 = sum(ThetaHist_c2p1, 2);
% PhiHist_c2p1 = sum(PhiHist_c2p1, 2);
% ThetaHist_c2p1 = ThetaHist_c2p1/sum(ThetaHist_c2p1);
% PhiHist_c2p1 = PhiHist_c2p1/sum(PhiHist_c2p1);
% BoltzDist_c2p1 = exp(c2p1*sqrt(3/8)*sin(2*ThetaBins));
% BoltzInt_c2p1 = sum(BoltzDist_c2p1.*sin(ThetaBins));
% BoltzDist_c2p1 = BoltzDist_c2p1.*sin(ThetaBins)./BoltzInt_c2p1;
% rmsd_c2p1 = sqrt(sum((ThetaHist_c2p1 - BoltzDist_c2p1).^2)/nBins);

Hist2D_c2p2 = sum(Hist2D_c2p2, 3);
ThetaHist_c2p2 = sum(ThetaHist_c2p2, 2);
PhiHist_c2p2 = sum(PhiHist_c2p2, 2);
ThetaHist_c2p2 = ThetaHist_c2p2/sum(ThetaHist_c2p2);
PhiHist_c2p2 = PhiHist_c2p2/sum(PhiHist_c2p2);
BoltzDist_c2p2 = exp(c2p2*sqrt(3/8)*sin(ThetaBins).^2);
BoltzInt_c2p2 = sum(BoltzDist_c2p2.*sin(ThetaBins));
BoltzDist_c2p2 = BoltzDist_c2p2.*sin(ThetaBins)./BoltzInt_c2p2;
rmsd_c2p2 = sqrt(sum((ThetaHist_c2p2 - BoltzDist_c2p2).^2)/nBins);

% ChiSquare = sum(((ThetaHist - BoltzDist).^2)./ThetaHist);

% This seems like a loose condition and should be investigated further
if rmsd_c2p2 > 5e-3 ...
%    || rmsd_c2p1 > 1e-2 %...
%    || rmsd_c2p2 > 1e-2
  err = 1;
%   plot(ThetaBins, ThetaHist_c20, ThetaBins, BoltzDist_c20)
%   plot(ThetaBins, ThetaHist_c2p1, ThetaBins, BoltzDist_c2p1)
%   plot(ThetaBins, ThetaHist_c2p2, ThetaBins, BoltzDist_c2p2)
%     bar3(ThetaBins(1:end-1),  Hist2D_c2p1)
%     bar3(ThetaBins(1:end-1),  Hist2D_c2p2)
else  
  err = 0;
%   plot(ThetaBins, ThetaHist_c20, ThetaBins, BoltzDist_c20)
%   plot(ThetaBins, ThetaHist_c2p1, ThetaBins, BoltzDist_c2p1)
%   plot(ThetaBins, ThetaHist_c2p2, ThetaBins, BoltzDist_c2p2)
%     bar3(ThetaBins(1:end-1),  Hist2D_c2p2)
end

data = [];

end
