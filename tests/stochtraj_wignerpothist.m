function ok = test()

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

nBins = 50;

% L,M,K,lambda to test
LMK = [1,1,0; 1,0,1; 1,1,1; 2,0,0; 2,0,1; 2,0,2];
lambda = 4; % same lambda for all LMK

alphaBins = linspace(0, 2*pi, nBins+1);
betaBins = linspace(0, pi, nBins/2+1);
gammaBins = linspace(0, 2*pi, nBins+1);

% form grids to be acted upon by Boltzmann distribution function expressions
[Agrid,Bgrid,Ggrid] = ndgrid(alphaBins,betaBins,gammaBins);

N = round(nSteps/2);

for j = 1:size(LMK,1)
  Sys.Potential = [LMK(j,:) lambda];
  [~, ~, qTraj] = stochtraj_diffusion(Sys,Par);  % extract quaternions from trajectories
  
  % pre-allocate array for 3D histograms
  Hist3D = zeros(nBins+1, nBins/2+1, nBins+1, nTraj);
  
  for iTraj=1:nTraj
    % discard first half of each trajectory
    [alpha, beta, gamma] = quat2euler(qTraj(:,N:end,iTraj),'active');
    % calculate 3D histogram using function obtained from Mathworks File Exchange
    Hist3D(:,:,:,iTraj) = histcnd([alpha;beta;gamma].',...
                                          {alphaBins,betaBins,gammaBins});
  end
  
  Hist3D = mean(Hist3D, 4);  % average over all trajectories
  Hist3D = Hist3D/sum(Hist3D(:));  % normalize
  
  rmsd = calc_rmsd(Sys.Potential(1,4),Sys.Potential(1,1:3),Hist3D,...
                   Agrid,Bgrid,Ggrid);
  
  ok = rmsd<5e-3 && all(~isnan(Hist3D(:)));
  % numerical result does not match analytical result
  if ~ok, break; end
  
end

end

% Helper function to compare numerical result with analytic expression
% -------------------------------------------------------------------------

function rmsd = calc_rmsd(lambda, LMK, Hist3D, Agrid, Bgrid, Ggrid)

% convert LMK indices to a string
LMKstr = num2str(LMK(:));
LMKstr = strcat(LMKstr(1),LMKstr(2),LMKstr(3));

Re = real(lambda);
Im = imag(lambda);

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

end
