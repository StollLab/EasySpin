function [err,data] = test(opt,olddata)

% Total transition intensity should be independent of orientation of k
% vector relative to B0, no matter what the crystal orientation and the
% polarization angle are.

clear Sys Exp
Sys.S = 1;
Sys.D = -3*30e3*[1 0.01];

alpha = rand*pi;

Sys.S = 1;
Sys.D = 30e3*2*[1 0.1];
Exp.CrystalOrientation = rand(1,3)*pi;
Opt.Transitions = [1 2; 1 3; 2 3];
beta = linspace(0,pi/2);

for a = numel(beta):-1:1
  Exp.Mode = [beta(a) alpha];
  [dum,Intensity(:,a)] = resfreqs_matrix(Sys,Exp,Opt);
end

Intensity = abs(sum(Intensity));

err = min(Intensity) < 0.99999*max(Intensity);
data = [];
