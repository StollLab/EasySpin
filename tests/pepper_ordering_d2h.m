function ok = test(opt)

% Spin center with D2h hamiltonian symmetry
Sys.g = [2 2.1 2.2];
Sys.lw = 1;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.Range = [300 345];

Opt.GridSize = 10;
Opt.GridSymmetry = 'D2h';

% width of Gaussian along beta and gamma
wbeta = deg2rad(8);
wgamma = deg2rad(8);

% center of Gaussian in first octant
beta0 = -deg2rad(45);
gamma0 = -deg2rad(45);

% center moved to other octants by D2h symmetry operations
beta0_ =  [beta0  beta0     beta0     beta0       pi-beta0  pi-beta0  pi-beta0  pi-beta0];
gamma0_ = [gamma0 pi-gamma0 pi+gamma0 2*pi-gamma0 gamma0    pi-gamma0 pi+gamma0 2*pi-gamma0]; 

for k = 1:numel(gamma0_)
  % use simple gaussian distribution function
  Exp.Ordering = @(beta,gamma) exp(-(gamma-gamma0_(k)).^2/wgamma^2-(beta-beta0_(k)).^2/wbeta^2);
  [B,spc(k,:)] = pepper(Sys,Exp,Opt);
end

for k = 1:numel(gamma0_)-1
  ok(k) = areequal(spc(k+1,:),spc(1,:),1e-10,'rel');
end

if opt.Display
  subplot(2,1,1)
  plot(B,spc);
  subplot(2,1,2)
  plot(B,spc-spc(1,:));
end
