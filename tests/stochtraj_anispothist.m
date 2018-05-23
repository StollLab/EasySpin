function [err,data] = test(opt,olddata)
% Check that using stochtraj with an anisotropic orienting potential 
% generates a proper distribution of orientations by 
% (1) calculating a histogram of Euler angles for each trajectory, and
% (2) averaging over all trajectories and comparing with the corresponding
%     Boltzmann distribution function, i.e. \exp(-U(\omega)/kT)

Sys.tcorr = 10e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(500*Sys.tcorr/Par.dt);
Par.nTraj = 400;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

nBins = 50;

% L,M,K indices to test
LMK = [1,1,0;
       1,0,1;
       1,1,1;
       2,0,0;
       2,0,1;
       2,0,2];
     
ReLambda = 2*rand(size(LMK,1),1)-1;  % real part of ordering coefficient lambda
ImLambda = 2*rand(size(LMK,1),1)-1;  % imaginary part

AlphaBins = linspace(-pi, pi, nBins);
BetaBins = linspace(0, pi, nBins/2);
GammaBins = linspace(-pi, pi, nBins);

% form grids to be acted upon by Boltzmann distribution function
% expressions
[Agrid,Bgrid,Ggrid] = meshgrid(AlphaBins,BetaBins,GammaBins);
     
err = 0;

N = round(nSteps/2);

for j=1:size(LMK,1)
  Sys.Coefs = [4, 0];
  Sys.LMK = LMK(j,:);
  [~, q] = stochtraj(Sys,Par);  % extract quaternions from trajectories
  
  % pre-allocate array for 3D histograms
  % note that the output will be of size (nBins,nBins,nBins), rather than the 
  % typical (nBins-1,nBins-1,nBins-1) size output from MATLAB's hist
  % functions
  Hist3D = zeros(nBins, nBins/2, nBins, nTraj);
  
  for iTraj=1:nTraj
    % use a "burn-in method" by taking last half of each trajectory
    [alpha, beta, gamma] = quat2euler(q(:,iTraj,N:end));
    alpha = squeeze(alpha);
    beta = squeeze(beta);
    gamma = squeeze(gamma);
    % calculate 3D histogram using function obtained from Mathworks File Exchange
    [Hist3D(:,:,:,iTraj),dummy] = histcnd([alpha,beta,gamma],...
                                          {AlphaBins,BetaBins,GammaBins});
  end

  Hist3D = mean(Hist3D, 4);  % average over all trajectories
  Hist3D = Hist3D/sum(Hist3D(:));  % normalize
  
%   slice(Agrid, Bgrid, Ggrid, Hist3D, 0, pi/2, 0)
%   xlabel('alpha')
%   ylabel('beta')
%   zlabel('gamma')
%   colormap hsv
  
  Hist3D = permute(Hist3D, [2, 1, 3]);  % first two dims are incorrectly 
                                        % ordered by histcnd
  
%   slice(Agrid, Bgrid, Ggrid, Hist3D, 0, pi/2, 0)
%   xlabel('alpha')
%   ylabel('beta')
%   zlabel('gamma')
%   colormap hsv
  
  rmsd = calc_rmsd(Sys.Coefs,Sys.LMK,Hist3D,Agrid,Bgrid,Ggrid);
  
  if rmsd>5e-3||any(isnan(Hist3D(:)))
    % numerical result does not match analytical result
    err = 1;
    break
  end

end

data = [];


% Helper function to compare numerical result with analytic expression
% -------------------------------------------------------------------------

function rmsd = calc_rmsd(Coefs, LMK, Hist3D, Agrid, Bgrid, Ggrid)

% convert LMK indices to a string
LMKstr = num2str(LMK(:));
LMKstr = strcat(LMKstr(1),LMKstr(2),LMKstr(3));

Re = Coefs(1);
Im = Coefs(2);

switch LMKstr
  case '110'
    BoltzDist = exp(1/sqrt(2)*sin(Bgrid).*(Re*cos(Ggrid)-Im*sin(Ggrid)));
  case '101'
    BoltzDist = exp(1/sqrt(2)*sin(Bgrid).*(Re*cos(Agrid)+Im*sin(Agrid)));
  case '111'
    BoltzDist = exp(-cos(Bgrid/2).^2.*(Re*cos(Agrid+Ggrid)-Im*sin(Agrid+Ggrid)));
  case '200'
    BoltzDist = exp(Re*1/2*(3*cos(Bgrid).^2 - 1));
  case '201'
    BoltzDist = exp(sqrt(3/8)*sin(2*Bgrid).*(Im*sin(Agrid)-Re*cos(Agrid)));
  case '202'
    BoltzDist = exp(sqrt(3/8)*sin(Bgrid).^2.*(Re*cos(2*Agrid)-Im*sin(2*Agrid)));
end

BoltzInt = sum(reshape(BoltzDist.*sin(Bgrid),1,numel(Bgrid)));
BoltzDist = BoltzDist.*sin(Bgrid)/BoltzInt;

% figure(1)
% slice(Agrid, Bgrid, Ggrid, Hist3D, 0, pi/2, 0)
% xlabel('alpha')
% ylabel('beta')
% zlabel('gamma')
% colormap hsv
% 
% figure(2)
% slice(Agrid, Bgrid, Ggrid, BoltzDist, 0, pi/2, 0)
% xlabel('alpha')
% ylabel('beta')
% zlabel('gamma')
% colormap hsv
% colorbar

residuals = Hist3D - BoltzDist;
rmsd = sqrt(mean(residuals(:).^2));

end

end
