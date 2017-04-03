function [err,data] = test(opt,olddata)
% Check that using stochtraj with an anisotropic orienting potential 
% generates a proper distribution of orientations by 
% (1) calculating a histogram of Euler angles for each trajectory, and
% (2) averaging over all trajectories and comparing with the corresponding
%     Boltzmann distribution function, i.e. \exp(-U(\omega)/kT)

Sys.tcorr = 10*rand()*1e-9;
Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
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

AlphaBins = linspace(-pi, pi, nBins)';
BetaBins = linspace(0, pi, nBins)';
GammaBins = linspace(-pi, pi, nBins)';

% form grids to be acted upon by Boltzmann distribution function
% expressions
[Agrid,Bgrid,Ggrid] = meshgrid(AlphaBins,BetaBins,GammaBins);

% pre-allocate array for 3D histograms
% note that the output will be of size (nBins,nBins,nBins), rather than the 
% typical (nBins-1,nBins-1,nBins-1) size output from MATLAB's hist
% functions
Hist3D = zeros(nBins, nBins, nBins, nTraj);
     
err = 0;

for j=1:size(LMK,1)
  Sys.Coefs = [ReLambda(j), ImLambda(j)];
  Sys.LMK = LMK(j,:);
  [~, ~, q] = stochtraj(Sys,Par);  % extract quaternions from trajectories
  
  for iTraj=1:nTraj
    % use a "burn-in method" by taking last half of each trajectory
    [alpha, beta, gamma] = quat2euler(q(:,iTraj,round(nSteps/2):end));
    % calculate 3D histogram using function obtained from Mathworks File Exchange
    [Hist3D(:,:,:,iTraj),~] = histcnd([alpha,beta,gamma],{AlphaBins',BetaBins',GammaBins'});
  end

  Hist3D = sum(Hist3D, 4);  % average over all trajectories
  Hist3D = Hist3D/sum(reshape(Hist3D,1,numel(Hist3D)));  % normalize
  
  if calc_rmsd(Sys.Coefs,Sys.LMK,Hist3D,Agrid,Bgrid,Ggrid)>5e-3
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

Re = Coefs(2);
Im = Coefs(2);

switch LMKstr
  case '110'
    BoltzDist = exp(1/sqrt(2)*sin(Bgrid).*(Re*cos(Ggrid)-Im*sin(Ggrid)));
  case '101'
    BoltzDist = exp(1/sqrt(2)*sin(Bgrid).*(Re*cos(Agrid)-Im*sin(Agrid)));
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
rmsd = sqrt(sum(reshape((Hist3D - BoltzDist).^2,1,numel(Hist3D))));

end

end
