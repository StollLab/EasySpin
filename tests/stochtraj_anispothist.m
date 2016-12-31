function [err,data] = test(opt,olddata)
% Check that using stochtraj with anisotropic diffusion generates a
% proper distribution of orientations

Sys.tcorr = 10*rand()*1e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;
% Par.beta = 0;%pi*(2*rand()-1);
% Par.alpha = 0;%2*pi*(2*rand()-1);
% Par.gamma = 0;%2*pi*(2*rand()-1);

nTraj = Par.nTraj;
nSteps = Par.nSteps;

nBins = 50;

c20 = 4;
Sys.Coefs = [c20, c20];
Sys.LMK = [2, 0, 0];
[~, ~, q_c20] = stochtraj(Sys,Par);

Rec2p1 = 4;
Imc2p1 = 4;
Sys.Coefs = [Rec2p1, Imc2p1];
Sys.LMK = [2, 0, 1];
[~, ~, q_c2p1] = stochtraj(Sys,Par);

Rec2p2 = 4;
Imc2p2 = 4;
Sys.Coefs = [Rec2p2, Imc2p2];
Sys.LMK = [2, 0, 2];
[~, ~, q_c2p2] = stochtraj(Sys,Par);


ThetaBins = linspace(0, pi, nBins)';
PhiBins = linspace(-pi, pi, nBins)';

Hist2D_c20 = zeros(nBins-1, nBins-1, nTraj);

Hist2D_c2p1 = zeros(nBins-1, nBins-1, nTraj);

Hist2D_c2p2 = zeros(nBins-1, nBins-1, nTraj);

for iTraj = 1:nTraj
  [phi, theta, ~] = quat2euler(q_c20(:,iTraj,:));
  [Hist2D_c20(:,:,iTraj),~,~] = histcounts2(theta,phi,ThetaBins,PhiBins);

  [phi, theta, ~] = quat2euler(q_c2p1(:,iTraj,:));
  [Hist2D_c2p1(:,:,iTraj),~,~] = histcounts2(theta,phi,ThetaBins,PhiBins);

  [phi, theta, ~] = quat2euler(q_c2p2(:,iTraj,:));
  [Hist2D_c2p2(:,:,iTraj),~,~] = histcounts2(theta,phi,ThetaBins,PhiBins);

end

% ThetaHist_c20 = sum(ThetaHist_c20,2);
% PhiHist_c20 = sum(PhiHist_c20,2);
% ThetaHist_c20 = ThetaHist_c20/sum(ThetaHist_c20);
% PhiHist_c20 = PhiHist_c20/sum(PhiHist_c20);
Hist2D_c20 = sum(Hist2D_c20, 3);
Hist2D_c20 = Hist2D_c20/sum(sum(Hist2D_c20,1),2);
[X,Y] = meshgrid((ThetaBins(2:end)+ThetaBins(1:end-1))/2, (PhiBins(2:end)+PhiBins(1:end-1))/2);
BoltzDist_c20 = exp(c20*1/2*(3*cos(X).^2 - 1));
BoltzInt_c20 = sum(sum(BoltzDist_c20.*sin(X),1),2);
BoltzDist_c20 = BoltzDist_c20.*sin(X)/BoltzInt_c20;
rmsd_c20 = sqrt(sum(sum((Hist2D_c20' - BoltzDist_c20).^2,1),2))

% ThetaHist_c2p1 = sum(ThetaHist_c2p1, 2);
% PhiHist_c2p1 = sum(PhiHist_c2p1, 2);
% ThetaHist_c2p1 = ThetaHist_c2p1/sum(ThetaHist_c2p1);
% PhiHist_c2p1 = PhiHist_c2p1/sum(PhiHist_c2p1);
Hist2D_c2p1 = sum(Hist2D_c2p1, 3);
Hist2D_c2p1 = Hist2D_c2p1/sum(sum(Hist2D_c2p1,1),2);
[X,Y] = meshgrid((ThetaBins(2:end)+ThetaBins(1:end-1))/2, (PhiBins(2:end)+PhiBins(1:end-1))/2);
BoltzDist_c2p1 = exp(sqrt(3/8)*sin(2*X).*(Imc2p1*sin(Y)-Rec2p1*cos(Y)));
BoltzInt_c2p1 = sum(sum(BoltzDist_c2p1.*sin(X),1),2);
BoltzDist_c2p1 = BoltzDist_c2p1.*sin(X)/BoltzInt_c2p1;
rmsd_c2p1 = sqrt(sum(sum((Hist2D_c2p1' - BoltzDist_c2p1).^2,1),2))

% ThetaHist_c2p2 = sum(ThetaHist_c2p2, 2);
% PhiHist_c2p2 = sum(PhiHist_c2p2, 2);
% ThetaHist_c2p2 = ThetaHist_c2p2/sum(ThetaHist_c2p2);
% PhiHist_c2p2 = PhiHist_c2p2/sum(PhiHist_c2p2);
Hist2D_c2p2 = sum(Hist2D_c2p2, 3);
Hist2D_c2p2 = Hist2D_c2p2/sum(sum(Hist2D_c2p2,1),2);
[X,Y] = meshgrid((ThetaBins(2:end)+ThetaBins(1:end-1))/2, (PhiBins(2:end)+PhiBins(1:end-1))/2);
BoltzDist_c2p2 = exp(sqrt(3/8)*sin(X).^2.*(Rec2p2*cos(2*Y)-Imc2p2*sin(2*Y)));
BoltzInt_c2p2 = sum(sum(BoltzDist_c2p2.*sin(X),1),2);
BoltzDist_c2p2 = BoltzDist_c2p2.*sin(X)/BoltzInt_c2p2;
rmsd_c2p2 = sqrt(sum(sum((Hist2D_c2p2' - BoltzDist_c2p2).^2,1),2))

% ChiSquare = sum(((ThetaHist - BoltzDist).^2)./ThetaHist);

% This seems like a loose condition and should be investigated further
if rmsd_c20 > 5e-3 ...
   || rmsd_c2p1 > 5e-3 ...
   || rmsd_c2p2 > 5e-3
  err = 1;
%   plot(ThetaBins, ThetaHist_c20, ThetaBins, BoltzDist_c20)
%   plot(ThetaBins, ThetaHist_c2p1, ThetaBins, BoltzDist_c2p1)
%   legend('Numeric','Analytic')
%   plot(ThetaBins, ThetaHist_c2p2, ThetaBins, BoltzDist_c2p2)
%   mesh(X, Y, Hist2D_c2p1' - BoltzDist_c2p1)
%   mesh(X, Y, Hist2D_c2p2' - BoltzDist_c2p2)
%   figure
%   subplot(2,1,1)  
%   mesh(X, Y, Hist2D_c2p2')
%   subplot(2,1,2)
%   mesh(X, Y, BoltzDist_c2p2)
%     bar3(ThetaBins(1:end-1),  Hist2D_c2p1)
%     bar3(ThetaBins(1:end-1),  Hist2D_c2p2)
else  
  err = 0;
end

data = [];

end
