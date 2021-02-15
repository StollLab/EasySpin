function ok = test(opt)

%=======================================================================
% Partially oriented system with axial spin parameters
%=======================================================================
Sys = struct('g',[2, 2, 2.2],'lw',1);
Exp = struct('mwFreq',9.5,'Range',[300 350]);
Opt = struct('Verbosity',0);

% Isotropic distribution
Exp.Ordering = [];
[B,spc0] = pepper(Sys,Exp,Opt);

% Isotropic distribution, user supplied
%Exp.Ordering = @(phi,theta) ones(size(phi));
Exp.Ordering = @ord_ones;
[B,spc1] = pepper(Sys,Exp,Opt);

% Gaussian distribution
%Exp.Ordering = @(phi,theta) gaussian(theta,0,pi/8);
Exp.Ordering = @ord_gaussian;
[B,spc2] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(B,spc0,B,spc1,B,spc2);
  legend('no ordering function','isotropic ordering function','Gaussian');
  legend boxoff
end

ok = areequal(spc0,spc1,1e-10,'rel');

end

function w = ord_ones(phi,theta)
w = ones(size(phi));
end

function w = ord_gaussian(phi,theta)
w = gaussian(theta,0,pi/8).*ones(size(phi));
end
