function ok = test(opt)

Sys.g = [2 2.1 2.2];
Sys.lw = 1;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.Range = [300 345];

Opt.GridSize = 10;
Opt.GridSymmetry = 'D2h';

% width of Gaussian along phi and theta
wphi = 8*pi/180;
wtheta = 8*pi/180;

% center of Gaussian in first octant
phi0 = 45*pi/180;
the0 = 45*pi/180;

% center moved to other octants by D2h symmetry operations
phi0_ = [phi0 pi-phi0 pi+phi0 2*pi-phi0 phi0    pi-phi0 pi+phi0 2*pi-phi0]; 
the0_ = [the0 the0    the0    the0      pi-the0 pi-the0 pi-the0 pi-the0];

for k = 1:numel(phi0_)
  % use simple gaussian distribution function
  Exp.Ordering = @(phi,theta)exp(-(phi-phi0_(k)).^2/wphi^2-(theta-the0_(k)).^2/wtheta^2);
  [B,spc(k,:)] = pepper(Sys,Exp,Opt);
end

for k = 1:numel(phi0_)-1
  ok(k) = areequal(spc(k+1,:),spc(1,:),1e-10,'rel');
end

if opt.Display
  subplot(2,1,1)
  plot(B,spc);
  subplot(2,1,2)
  plot(B,spc-spc(1,:));
end
