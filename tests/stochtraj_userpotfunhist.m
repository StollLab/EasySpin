function [err,data] = test(opt,olddata)
% Check that supplying a pseudopotential energy function to stochtraj with 
% generates a proper distribution of orientations

% generate Euler angle grids
N = 100;
agrid = linspace(-180, 180, N)/180*pi;
bgrid = linspace(-90, 90, N)/180*pi;
ggrid = linspace(-180, 180, N)/180*pi;
% [Agrid, Bgrid, Ggrid] = meshgrid(agrid, bgrid, ggrid);

fwhma = 60/pi*180;
fwhmb = 60/pi*180;
fwhmg = 60/pi*180;

% PotFun = gaussian(Agrid, 0, fwhma).*gaussian(Bgrid, 0, fwhmb) ...
%                                   .*gaussian(Ggrid, 0, fwhmg);
[funa, funb, func] = meshgrid(wrappedgaussian(agrid, 0, fwhma), ...
                              wrappedgaussian(bgrid, 0, fwhmb), ...
                              wrappedgaussian(ggrid, 0, fwhmg));
PotFun = funa.*funb.*func;

Sys.tcorr = 10*1e-9;
Sys.PseudoPotFun = PotFun;

Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

nBins = 100;

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


[~, ~, q] = stochtraj(Sys,Par);  % extract quaternions from trajectories

for iTraj=1:nTraj
  % use a "burn-in method" by taking last half of each trajectory
  [alpha, beta, gamma] = quat2euler(q(:,iTraj,round(nSteps/2):end));
  % calculate 3D histogram using function obtained from Mathworks File Exchange
  [Hist3D(:,:,:,iTraj),~] = histcnd([alpha,beta,gamma],{AlphaBins',BetaBins',GammaBins'});
end

Hist3D = sum(Hist3D, 4);  % average over all trajectories
Hist3D = Hist3D/sum(reshape(Hist3D,1,numel(Hist3D)));  % normalize
rmsd = calc_rmsd(PotFun,Hist3D,Agrid,Bgrid,Ggrid);

if rmsd>5e-10
  % numerical result does not match analytical result
  err = 1;
end

data = [];


% Helper function to compare numerical result with analytic expression
% -------------------------------------------------------------------------

function rmsd = calc_rmsd(PotFun, Hist3D, Agrid, Bgrid, Ggrid)

rmsd = sqrt(sum(reshape((Hist3D - PotFun).^2,[1,numel(Hist3D)])));

end

end
