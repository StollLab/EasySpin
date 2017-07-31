function [err,data] = test(opt,olddata)
% Check that supplying a pseudopotential energy function to stochtraj 
% generates a proper distribution of orientations

% generate Euler angle grids

abins = 40;
bbins = abins/2;
gbins = abins;

agrid = linspace(-180, 180, abins)/180*pi;
bgrid = linspace(0, 180, bbins)/180*pi;
ggrid = linspace(-180, 180, gbins)/180*pi;
% [Agrid, Bgrid, Ggrid] = meshgrid(agrid, bgrid, ggrid);

delta = 30;

fwhma = delta/180*pi;
fwhmb = delta/180*pi;
fwhmg = delta/180*pi;

% [funa, funb, fung] = meshgrid(gaussian(agrid, 0, fwhma), ...
%                               gaussian(bgrid, 0, fwhma), ...
%                               gaussian(ggrid, 0, fwhma));

[funa, funb, fung] = ndgrid(wrappedgaussian(agrid, 0, fwhma), ...
                              wrappedgaussian(bgrid, 0, fwhmb, [0,pi]), ...
                              wrappedgaussian(ggrid, 0, fwhmg));

PotFun = funa.*funb.*fung;

Sys.tcorr = 1e-9;
Sys.PseudoPotFun = PotFun;

Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(500*Sys.tcorr/Par.dt);
Par.nTraj = 400;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

AlphaBins = linspace(-pi, pi, abins);
BetaBins = linspace(0, pi, bbins);
GammaBins = linspace(-pi, pi, gbins);

% form grids to be acted upon by Boltzmann distribution function
% expressions
[Agrid,Bgrid,Ggrid] = ndgrid(AlphaBins,acos(linspace(1,-1,bbins)),GammaBins);

% pre-allocate array for 3D histograms
% note that the output will be of size (nBins,nBins,nBins), rather than the 
% typical (nBins-1,nBins-1,nBins-1) size output from MATLAB's hist
% functions
% Hist3D = zeros(abins, bbins, gbins, nTraj);
     
err = 0;

Par.Omega = [0;0;0];

[t, RTraj, qTraj] = stochtraj(Sys,Par);  % extract quaternions from trajectories

M = 500;

% hold on
% [X,Y,Z] = sphere(25);
% a = 0.99;
% h = surf(a*X, a*Y, a*Z);
% set(h, 'edgecolor','none');
% colormap(1,[0,0,1]);
% alpha 0.1
% 
% VecTraj = squeeze(RTraj(:,3,1,:));
% plot3(squeeze(VecTraj(1,1:M)), ...
%       squeeze(VecTraj(2,1:M)), ...
%       squeeze(VecTraj(3,1:M)), ...
%       'LineWidth', 1.5,...
%       'Color', 'red')
%     
% camlight left
% set(gca,'CameraViewAngle',6);
% view([45,45])
% daspect([1,1,1]);
% axis off
% hold off

% N = round(nSteps/2);
N = 1;

for iTraj=1:nTraj
  % use a "burn-in method" by taking last half of each trajectory
  [alph, beta, gamma] = quat2euler(qTraj(:,iTraj,N:end));
  alph = squeeze(alph);
  beta = squeeze(beta);
  gamma = squeeze(gamma);
  
  % calculate 3D histogram using function obtained from Mathworks File Exchange
  [Hist3D(:,:,:,iTraj),~] = histcnd([alph,beta,gamma],...
                                    {AlphaBins,BetaBins,GammaBins});
end

Hist3D = mean(Hist3D, 4);  % average over all trajectories
Hist3D = Hist3D/sum(Hist3D(:));  % normalize

% Hist3D = permute(Hist3D, [2, 1, 3]);  % first two dims are incorrectly 
%                                       % ordered by histcnd

idx = [2, 1, 3];

% figure(1)
% slice(permute(pi/2-Agrid, idx), ...
%       permute(Bgrid, idx), ...
%       permute(Ggrid, idx), ...
%       permute(Hist3D, idx), ...
%       0, pi/2, 0)
% xlabel('alpha')
% ylabel('beta')
% zlabel('gamma')
% colormap hsv
% 
% figure(2)
% slice(permute(Agrid, idx), ...
%       permute(Bgrid, idx), ...
%       permute(Ggrid, idx), ...
%       permute(PotFun, idx), 0, pi/2, 0)
% xlabel('alpha')
% ylabel('beta')
% zlabel('gamma')
% colormap hsv

% figure(2)
% plot3(unwrap(alpha), beta, unwrap(gamma))

rmsd = calc_rmsd(PotFun, Hist3D);

if rmsd>1e-2||any(isnan(Hist3D(:)))
  % numerical result does not match analytical result
  err = 1;
end

data = [];


% Helper function to compare numerical result with analytic expression
% -------------------------------------------------------------------------

function rmsd = calc_rmsd(PotFun, Hist3D)

residuals = Hist3D - PotFun;
rmsd = sqrt(mean(reshape(residuals.^2,[1,numel(Hist3D)])));

end

end
