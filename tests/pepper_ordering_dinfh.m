function ok = test(opt)

Sys.g = [2 2 2.2];
Sys.lw = 1;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.Range = [300 345];

Opt.GridSize = 10;
Opt.GridSymmetry = 'Dinfh';

% width of Gaussian along phi and theta
wbeta = 8*pi/180;
wgamma = 8*pi/180;

% center of Gaussian
beta0 = 45*pi/180;
gamma0 = 45*pi/180;

% center moved to other orientations by D2h symmetry operations
beta0_ = [beta0 beta0         pi-beta0  pi-beta0];
gamma0_ = [gamma0 gamma0+pi+pi/7    gamma0     gamma0+6*pi/5]; 

for k = 1:numel(gamma0_)
  % use simple gaussian distribution function
  Exp.Ordering = @(beta,gamma)exp(-(gamma-gamma0_(k)).^2/wgamma^2-(beta-beta0_(k)).^2/wbeta^2);
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
