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
Sys.ProbDensFun = pdf;

Par.dt = Sys.tcorr/20;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 400;

nTraj = Par.nTraj;
nSteps = Par.nSteps;

AlphaBins = linspace(-pi, pi, abins);
BetaBins = linspace(0, pi, bbins);
GammaBins = linspace(-pi, pi, gbins);

tcorr = Sys.tcorr;
nTraj = Par.nTraj;
nSteps = Par.nSteps;

[t, q] = stochtraj(Sys,Par);
R = quat2rotmat(q);

VecTraj = squeeze(R(:,3,:,:));

AutoCorrFFT = runprivate('autocorrfft', VecTraj.^2, 3, 1, 1);

AutoCorrFFT = squeeze(mean(AutoCorrFFT, 2));

N = round(nSteps/2);

AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(N:end));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-(1/tcorr)*t);

[k, c, yFit] = exponfit(t(1:N), AutoCorrFFT(1:N));
tauR = 1/k;

residuals = AutoCorrFFT(1:N) - yFit;
rmsd = sqrt(mean(residuals.^2));

if rmsd > 1e-2 || isnan(rmsd) || tcorr-tauR < 0
  err = 1;
  plot(t(1:N)/1e-6, AutoCorrFFT(1:N), t(1:N)/1e-6, analytic(1:N))
  xlabel('t (\mu s)')
else
  err = 0;
end

data = [];

end
