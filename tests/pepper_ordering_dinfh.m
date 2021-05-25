function ok = test(opt)

Sys.g = [2 2 2.2];
Sys.lw = 1;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.Range = [300 345];

Opt.GridSize = 10;
Opt.GridSymmetry = 'Dinfh';

% width of Gaussian along phi and theta
wphi = 8*pi/180;
wtheta = 8*pi/180;

% center of Gaussian
phi0 = 45*pi/180;
the0 = 45*pi/180;

% center moved to other orientations by D2h symmetry operations
phi0_ = [phi0 phi0+pi+pi/7    phi0     phi0+6*pi/5]; 
the0_ = [the0 the0         pi-the0  pi-the0];

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
