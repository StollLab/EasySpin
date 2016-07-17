% Signal averaging using the EWRLS adaptive filter method
%==========================================================================
% This example illustrates ewrls(), a function that computes fitered averages
% of multiple scans that have in general a better signal-to-noise ratio than
% the straighforward average. The function uses an exponentially weighted
% recursive least squares adaptive filter. The quality of the filtering
% depends on two parameters, the filter length p and the memory time lambda.
% Care has to be taken to pick values that do not distort the signal.

clear, clc, clf

% Generate noisy example data
%--------------------------------------------------------
nScans = 40;
nPoints = 1001;

t = linspace(-1,1,nPoints);

signal = gaussian(t,0.2,0.3,1) + gaussian(t,-0.2,0.15,1)/2;
signal = signal(:)/max(signal);
for k = 1:nScans
  datamatrix(:,k) = addnoise(signal,2,'n');
end

% run EWRLS filtering
%---------------------------------------------------------
p = 20; % filter length
lambda = 0.96; % memory factor
nPreAverage = 0;
delta = 100; % regularization parameter
ewrls(datamatrix,p,lambda,nPreAverage,delta,'fb');
