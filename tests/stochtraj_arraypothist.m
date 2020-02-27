function ok = test(opt)

% Check that supplying a pseudopotential energy function to stochtraj_diffusion 
% generates a proper distribution of orientations.

rng(1)

% generate Euler angle grids
N = 70;

na = N+1;
nb = N/2+1;
ng = N+1;

alphaGrid = linspace(0, 2*pi, na);
betaGrid = linspace(0, pi, nb);
gammaGrid = linspace(0, 2*pi, ng);

delta = 80;

fwhma = delta/180*pi;
fwhmb = delta/180*pi;
fwhmg = delta/180*pi;

[pdfa, pdfb, pdfg] = ...
  ndgrid(runprivate('wrappedgaussian', alphaGrid, pi/2, fwhma, [0,2*pi]), ...
         runprivate('wrappedgaussian', betaGrid, 0, fwhmb, [0,pi]), ...
         runprivate('wrappedgaussian', gammaGrid, 0, fwhmg, [0,2*pi]));

pdf1 = pdfa.*pdfb.*pdfg;

[pdfa, pdfb, pdfg] = ndgrid(runprivate('wrappedgaussian', alphaGrid, pi/4, fwhma, [0,2*pi]), ...
                            runprivate('wrappedgaussian', betaGrid, pi/2, fwhmb, [0,pi]), ...
                            runprivate('wrappedgaussian', gammaGrid, 3*pi/4, fwhmg, [0,2*pi]));

pdf2 = pdfa.*pdfb.*pdfg;
pdf = pdf1 + pdf2;
                          
pdf = pdf/sum(pdf(:));
PDF = pdf;
PDF(PDF<1e-10) = 1e-10;
U = -log(PDF);

Sys.tcorr = 10e-9;
Sys.Potential = U;

Par.dt = Sys.tcorr/20;
Par.nSteps = ceil(400*Sys.tcorr/Par.dt);
Par.nTraj = 400;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

% pre-allocate array for 3D histograms
Hist3D = zeros(na,nb,ng,nTraj);

% Generate quaternion trajectories
[~,~,qTraj] = stochtraj_diffusion(Sys,Par);
reverse = repmat(qTraj(1,:,:)<0,[4,1,1]);
qTraj(reverse) = -qTraj(reverse);

alphaBins = linspace(0, 2*pi, na);
betaBins = linspace(0, pi, nb);
gammaBins = linspace(0, 2*pi, ng);

N = round(nSteps/2);
for iTraj = 1:nTraj
  % discard first half of each trajectory
  [alpha, beta, gamma] = quat2euler(qTraj(:,N:end,iTraj),'active');
  
  % calculate 3D histogram using function obtained from Mathworks File Exchange
  [Hist3D(:,:,:,iTraj),~] = histcnd([alpha;beta;gamma].',...
                                    {alphaBins,betaBins,gammaBins});
end

Hist3D = mean(Hist3D,4);  % average over all trajectories
Hist3D = Hist3D/sum(Hist3D(:));  % normalize

pdf = exp(-U);

pdf = pdf.*sin(betaBins)/sum(pdf(:));

rmsd = calc_rmsd(pdf, Hist3D);

ok = rmsd<0.3 && all(~isnan(Hist3D(:)));
% error if numerical result does not match analytical result

if opt.Display
  subplot(1,2,1)
  slice(alphaBins, ...
        betaBins, ...
        gammaBins, ...
        permute(Hist3D, [2 1 3]), ...
        0, pi/2, 0)
  xlabel('alpha')
  ylabel('beta')
  zlabel('gamma')
  title('Histogram')
  colormap hsv

  subplot(1,2,2)
  slice(alphaBins, ...
        betaBins, ...
        gammaBins, ...
        permute(pdf, [2 1 3]), ...
        0, pi/2, 0)
  xlabel('alpha')
  ylabel('beta')
  zlabel('gamma')
  title('PDF from potential')
  colormap hsv
end


% Helper function to compare numerical result with analytic expression
% -------------------------------------------------------------------------

function rmsd = calc_rmsd(PotFun, Hist3D)

residuals = Hist3D - PotFun;
rmsd = sqrt(mean(residuals(:).^2))/max(PotFun(:));

end

end
