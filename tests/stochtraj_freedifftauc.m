function [err,data] = test(opt,olddata)

Par.tcorr = 10*rand()*1e-9;
Par.dt = 1.0e-9;
Par.nSteps = 2000;
Par.nTraj = 800;

tcorr = Par.tcorr;
nTraj = Par.nTraj;
nSteps = Par.nSteps;

[t, R] = stochtraj(Par);

VecTraj = squeeze(R(:, 3, :, :));

AutoCorr = zeros(nSteps, nTraj);

for iTraj = 1:nTraj
  AutoCorr(:, iTraj) = autocorrfft(VecTraj(3, :, iTraj).^2);
end

AutoCorr = sum(AutoCorr, 2)'/nTraj;


analytic = exp(-(1/tcorr)*t);

ChiSquare = sum(((AutoCorr - analytic).^2)./AutoCorr);

if ChiSquare > 1e-4
  err = 1;
else  
  err = 0;
end

data = [];

end
