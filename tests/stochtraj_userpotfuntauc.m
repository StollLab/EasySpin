function [err,data] = test(opt,olddata)
% Check that supplying a pseudopotential energy function to stochtraj  
% generates a proper rotational correlation time

% generate Euler angle grids
N = 100;
agrid = linspace(-180, 180, N)/180*pi;
bgrid = linspace(0, 180, N)/180*pi;
ggrid = linspace(-180, 180, N)/180*pi;
% [Agrid, Bgrid, Ggrid] = meshgrid(agrid, bgrid, ggrid);

fwhma = 60/pi*180;
fwhmb = 60/pi*180;
fwhmg = 60/pi*180;

% PotFun = gaussian(Agrid, 0, fwhma).*gaussian(Bgrid, 0, fwhmb) ...
%                                   .*gaussian(Ggrid, 0, fwhmg);
[funa, funb, func] = meshgrid(wrappedgaussian(agrid, 0, fwhma), ...
                              wrappedgaussian(bgrid, 0, fwhmb, [0,pi]), ...
                              wrappedgaussian(ggrid, 0, fwhmg));
PotFun = funa.*funb.*func;

Sys.tcorr = 10*1e-9;
Sys.PseudoPotFun = PotFun;

Par.dt = Sys.tcorr/10;
Par.nSteps = ceil(200*Sys.tcorr/Par.dt);
Par.nTraj = 200;
Par.Omega = [  pi*(2*rand()-1); 
             2*pi*(2*rand()-1);
             2*pi*(2*rand()-1) ];

tcorr = Sys.tcorr;
nTraj = Par.nTraj;
nSteps = Par.nSteps;

[t, R] = stochtraj(Sys,Par);

VecTraj = squeeze(R(:,3,:,:));

AutoCorrFFT = zeros(nTraj, nSteps);

for iTraj = 1:nTraj
  AutoCorrFFT(iTraj,:) = autocorrfft(squeeze(VecTraj(:,iTraj,:).^2), 1);
end

N = round(nSteps/2);

AutoCorrFFT = mean(AutoCorrFFT, 1).';
AutoCorrFFT = AutoCorrFFT-mean(AutoCorrFFT(N:end));
AutoCorrFFT = AutoCorrFFT/max(AutoCorrFFT);

analytic = exp(-(1/tcorr)*t);

% ChiSquare = sum(((AutoCorrFFT - analytic).^2)./AutoCorrFFT)
residuals = AutoCorrFFT - analytic;
rmsd = sqrt(mean(residuals(1:N).^2));

if rmsd > 1e-2 || isnan(rmsd)
  err = 1;
  plot(t(1:N)/1e-6, AutoCorrFFT(1:N), t(1:N)/1e-6, analytic(1:N))
  xlabel('t (\mu s)')
else
  err = 0;
end

data = [];

end
