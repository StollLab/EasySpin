function [err,data] = test(opt,olddata)
% Check that supplying a pseudopotential energy function to stochtraj  
% generates a proper rotational correlation time

idx = [2, 1, 3];  % for permuting dimensions to go between ngrid and 
                  % meshgrid ordering

% generate Euler angle grids
N = 100;

abins = N;
bbins = N;
gbins = N;

alphaGrid = linspace(0, 2*pi, abins);
betaGrid = linspace(0, pi, bbins);
gammaGrid = linspace(0, 2*pi, gbins);

delta = 80;

fwhma = delta/180*pi;
fwhmb = delta/180*pi;
fwhmg = delta/180*pi;

[pdfa, pdfb, pdfg] = ndgrid(runprivate('wrappedgaussian', alphaGrid, 0, fwhma, [0,2*pi]), ...
                            runprivate('wrappedgaussian', betaGrid, 0, fwhmb, [0,pi]), ...
                            runprivate('wrappedgaussian', gammaGrid, 0, fwhmg, [0,2*pi]));

pdf = pdfa.*pdfb.*pdfg;

pdf = pdf/sum(pdf(:));

Sys.tcorr = 10e-9;
Sys.Potential = -log(pdf);

Par.dt = Sys.tcorr/20;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;

tcorr = Sys.tcorr;
nSteps = Par.nSteps;

[t, RTraj] = stochtraj_diffusion(Sys,Par);

VecTraj = squeeze(RTraj(:,3,:,:));

AutoCorrFFT = runprivate('autocorrfft', VecTraj.^2, 2, 1, 1);

AutoCorrFFT = squeeze(mean(AutoCorrFFT, 3));

N = round(nSteps/2);

AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(N:end));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-t/tcorr);

[k, ~, yFit] = exponfit(t(1:N), AutoCorrFFT(1:N));
tauR = 1/k;

residuals = AutoCorrFFT(1:N) - yFit;
rmsd = sqrt(mean(residuals.^2));

err = rmsd > 1e-2 || isnan(rmsd) || tcorr-tauR < 0;
if opt.Display
  x = t(1:N)/1e-6;
  plot(x, AutoCorrFFT(1:N), x, analytic(1:N));
  legend('autcorrFFT','analytical');
  xlabel('t (\mus)')
end

data = [];

end
