function [err,data] = test(opt,olddata)
% Check that using stochtraj_diffusion with an anisotropic orienting 
% potential generates a proper distribution of orientations by 
% (1) calculating a histogram of Euler angles for each trajectory, and
% (2) averaging over all trajectories and comparing with the corresponding
%     Boltzmann distribution function, i.e. \exp(-U(\omega)/kT)

Sys.tcorr = 10e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(500*Sys.tcorr/Par.dt);
Par.nTraj = 400;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

pidx = [2,1,3];

nBins = 50;

% L,M,K indices to test
LMK = [1,1,0;
       1,0,1;
       1,1,1;
       2,0,0;
       2,0,1;
       2,0,2];
% LMK = [2,0,0;
%        2,0,2];
lambda = 4; % same lambda for all LMK

alphaBins = linspace(0, 2*pi, nBins+1);
betaBins = linspace(0, pi, nBins/2+1);
gammaBins = linspace(0, 2*pi, nBins+1);

% form grids to be acted upon by Boltzmann distribution function
% expressions
% [Agrid,Bgrid,Ggrid] = meshgrid(AlphaBins,BetaBins,GammaBins);
[Agrid,Bgrid,Ggrid] = ndgrid(alphaBins,betaBins,gammaBins);
     
err = 0;

N = round(nSteps/2);

for j=1:size(LMK,1)
  Sys.Potential = [LMK(j,:) lambda];
  [~, ~, qTraj] = stochtraj_diffusion(Sys,Par);  % extract quaternions from trajectories
  
  % pre-allocate array for 3D histograms
  % note that the output will be of size (nBins,nBins,nBins), rather than the 
  % typical (nBins-1,nBins-1,nBins-1) size output from MATLAB's hist
  % functions
  Hist3D = zeros(nBins+1, nBins/2+1, nBins+1, nTraj);
  
  for iTraj=1:nTraj
    % use a "burn-in method" by taking last half of each trajectory
%     [alpha, beta, gamma] = quat2euler(qTraj(:,iTraj,N:end));
    [alpha, beta, gamma] = quat2euler(qTraj(:,iTraj,N:end),'active');
    alpha = squeeze(alpha);
    beta = squeeze(beta);
    gamma = squeeze(gamma);
    % calculate 3D histogram using function obtained from Mathworks File Exchange
    [Hist3D(:,:,:,iTraj),dummy] = histcnd([alpha,beta,gamma],...
                                          {alphaBins,betaBins,gammaBins});
  end

  Hist3D = mean(Hist3D, 4);  % average over all trajectories
  Hist3D = Hist3D/sum(Hist3D(:));  % normalize
  
  rmsd = calc_rmsd(Sys.Potential(1,4),Sys.Potential(1,1:3),Hist3D,...
                   Agrid,Bgrid,Ggrid,alphaBins,betaBins,gammaBins);
                 
  
  if rmsd>5e-3||any(isnan(Hist3D(:)))
    % numerical result does not match analytical result
    err = 1;
    break
  end

end

data = [];


% Helper function to compare numerical result with analytic expression
% -------------------------------------------------------------------------

function rmsd = calc_rmsd(lambda, LMK, Hist3D, Agrid, Bgrid, Ggrid,alphaBins,betaBins,gammaBins)

% convert LMK indices to a string
LMKstr = num2str(LMK(:));
LMKstr = strcat(LMKstr(1),LMKstr(2),LMKstr(3));

Re = real(lambda);
Im = imag(lambda);

% BoltzDist = exp(2*real(lambda*conj(wignerd(LMK([1,3,2]),Agrid,Bgrid,Ggrid))));
% BoltzDist = exp(2*real(lambda*conj(wignerd(LMK,Agrid,Bgrid,Ggrid))));
% BoltzDist = exp(2*real(lambda*conj(permute(wignerd(LMK,Agrid,Bgrid,Ggrid),[3,2,1]))));

switch LMKstr
  case '110'
%     BoltzDist = exp(1/sqrt(2)*sin(Bgrid).*(Re*cos(Ggrid)-Im*sin(Ggrid)));
    U = sqrt(2)*sin(Bgrid).*(Re*cos(Agrid)+Im*sin(Agrid));
%     BoltzDist = exp(2*real(lambda*conj(wignerd([1,1,0],Agrid,Bgrid,Ggrid))));
  case '101'
%     BoltzDist = exp(1/sqrt(2)*sin(Bgrid).*(Re*cos(Agrid)+Im*sin(Agrid)));
    U = -sqrt(2)*sin(Bgrid).*(Re*cos(Ggrid)+Im*sin(Ggrid));
%     BoltzDist = exp(2*real(lambda*conj(wignerd([1,0,1],Agrid,Bgrid,Ggrid))));
  case '111'
%     BoltzDist = exp(-cos(Bgrid/2).^2.*(Re*cos(Agrid+Ggrid)-Im*sin(Agrid+Ggrid)));
    U = -2*cos(Bgrid/2).^2.*(Re*cos(Agrid+Ggrid)+Im*sin(Agrid+Ggrid));
  case '200'
    U = -Re*(3*cos(Bgrid).^2 - 1);
  case '201'
    U = -sqrt(6)*cos(Bgrid).*sin(Bgrid).*(Re*cos(Ggrid)+Im*sin(Ggrid));
  case '202'
    U = -sqrt(3/2)*sin(Bgrid).^2.*(Re*cos(2*Ggrid)+Im*sin(2*Ggrid));
end

BoltzDist = exp(-U);

BoltzInt = sum(reshape(BoltzDist.*sin(Bgrid),1,numel(Bgrid)));
BoltzDist = BoltzDist.*sin(Bgrid)/BoltzInt;

residuals = Hist3D - BoltzDist;
rmsd = sqrt(mean(residuals(:).^2));

% subplot(1,2,1)
%   slice(alphaBins, ...
%         betaBins, ...
%         gammaBins, ...
%         permute(Hist3D, [2,1,3]), ...
%         0, pi/2, 0)
%   xlabel('alpha')
%   ylabel('beta')
%   zlabel('gamma')
%   title('Histogram')
%   colormap hsv
% 
%   subplot(1,2,2)
%   slice(alphaBins, ...
%         betaBins, ...
%         gammaBins, ...
%         permute(BoltzDist, [2,1,3]), ...
%         0, pi/2, 0)
%   xlabel('alpha')
%   ylabel('beta')
%   zlabel('gamma')
%   title('PDF from potential')
%   colormap hsv

end

end
