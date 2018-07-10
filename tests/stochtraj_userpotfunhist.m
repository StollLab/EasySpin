function [err,data] = test(opt,olddata)
% Check that supplying a pseudopotential energy function to stochtraj_diffusion 
% generates a proper distribution of orientations

idx = [2, 1, 3];  % for permuting dimensions to go between ngrid and 
                  % meshgrid ordering

% generate Euler angle grids
N = 100;

abins = N;
bbins = N;
gbins = N;

alphaGrid = linspace(-pi, pi, abins);
betaGrid = linspace(0, pi, bbins);
gammaGrid = linspace(-pi, pi, gbins);

delta = 80;

fwhma = delta/180*pi;
fwhmb = delta/180*pi;
fwhmg = delta/180*pi;

[pdfa, pdfb, pdfg] = ndgrid(runprivate('wrappedgaussian', alphaGrid, 0, fwhma), ...
                            runprivate('wrappedgaussian', betaGrid, 0, fwhmb, [0,pi]), ...
                            runprivate('wrappedgaussian', gammaGrid, 0, fwhmg));

pdf = pdfa.*pdfb.*pdfg;

pdf = pdf/sum(pdf(:));

Sys.tcorr = 10e-9;
Sys.Potential = -log(pdf);

Par.dt = Sys.tcorr/20;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

AlphaBins = linspace(-pi, pi, abins);
BetaBins = linspace(0, pi, bbins);
GammaBins = linspace(-pi, pi, gbins);

% pre-allocate array for 3D histograms
% note that the output will be of size (nBins,nBins,nBins), rather than the 
% typical (nBins-1,nBins-1,nBins-1) size output from MATLAB's hist
% functions
Hist3D = zeros(abins, bbins, gbins, nTraj);
     
err = 0;

[t, ~, qTraj] = stochtraj_diffusion(Sys,Par);  % extract quaternions from trajectories

M = 500;

N = round(nSteps/2);

alphaBins = linspace(-pi, pi, gbins);
betaBins = pi-acos(linspace(-1, 1, bbins));
gammaBins = linspace(-pi, pi, gbins);

for iTraj=1:nTraj
  % use a "burn-in method" by taking last half of each trajectory
  [alpha, beta, gamma] = quat2euler(qTraj(:,iTraj,N:end),'active');
  alpha = squeeze(alpha);
  beta = squeeze(beta);
  gamma = squeeze(gamma);
  
  % calculate 3D histogram using function obtained from Mathworks File Exchange
  [Hist3D(:,:,:,iTraj),~] = histcnd([alpha,beta,gamma],...
                                    {alphaBins,betaBins,gammaBins});
end

Hist3D = mean(Hist3D, 4);  % average over all trajectories
Hist3D = Hist3D/sum(Hist3D(:));  % normalize

[pdfa, pdfb, pdfg] = ndgrid(runprivate('wrappedgaussian', alphaBins, 0, fwhma), ...
                            runprivate('wrappedgaussian', betaBins, 0, fwhmb, [0,pi]).*sin(betaBins), ...
                            runprivate('wrappedgaussian', gammaBins, 0, fwhmg));

pdf = pdfa.*pdfb.*pdfg;

pdf = pdf/sum(pdf(:));

rmsd = calc_rmsd(pdf, Hist3D);

if rmsd>2e-2||any(isnan(Hist3D(:)))
  % numerical result does not match analytical result
  err = 1;
  
  subplot(1,2,1)
  slice(alphaBins, ...
        betaBins, ...
        gammaBins, ...
        permute(Hist3D, idx), ...
        0, pi/2, 0)
  xlabel('alpha')
  ylabel('beta')
  zlabel('gamma')
  colormap hsv

  subplot(1,2,2)
  slice(alphaBins, ...
        betaBins, ...
        gammaBins, ...
        permute(pdf, idx), ...
        0, pi/2, 0)
  xlabel('alpha')
  ylabel('beta')
  zlabel('gamma')
  colormap hsv
end

data = [];


% Helper function to compare numerical result with analytic expression
% -------------------------------------------------------------------------

function rmsd = calc_rmsd(PotFun, Hist3D)

residuals = Hist3D - PotFun;
rmsd = mean(abs(residuals(:)))/max(abs(residuals(:)));

end

end
