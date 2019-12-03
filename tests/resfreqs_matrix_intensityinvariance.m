function err = test(opt)

% Total transition intensity should be independent of orientation of k
% vector relative to B0, no matter what the crystal orientation and the
% polarization angle are.

Sys.S = 1;
Sys.D = -3*30e3*[1 0.01];

alpha = 13*pi/37;

Sys.S = 1;
Sys.D = 30e3*2*[1 0.1];
Exp.CrystalOrientation = [0.71 0.475 0.11]*pi;
Opt.Transitions = [1 2; 1 3; 2 3];
beta = linspace(0,pi/2,31);

for a = numel(beta):-1:1
  Exp.Mode = [beta(a) alpha];
  [~,Intensity(:,a)] = resfreqs_matrix(Sys,Exp,Opt);
end

Intensity = abs(sum(Intensity));

err = min(Intensity)/max(Intensity) < 0.99999;
