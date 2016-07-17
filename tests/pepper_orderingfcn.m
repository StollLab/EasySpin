function [err,data] = test(opt,olddata)

%=======================================================================
% Partially oriented system with axial spin parameters
%=======================================================================
Sys = struct('g',[2, 2, 2.2],'lw',1);
Exp = struct('mwFreq',9.5,'Range',[300 350]);
Opt = struct('Verbosity',0);

% Isotropic distribution
Exp.Ordering = [];
[x,spc0] = pepper(Sys,Exp,Opt);

% Isotropic distribution, user supplied
%Exp.Ordering = @(phi,theta) ones(size(phi));
Exp.Ordering = @ord_ones;
[x,spc1] = pepper(Sys,Exp,Opt);

% Gaussian distribution
%Exp.Ordering = @(phi,theta) gaussian(theta,0,pi/8);
Exp.Ordering = @ord_gaussian;
[x,spc2] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(x,spc0,x,spc1,'r');
  legend('no ordering function','isotropic ordering function');
  legend boxoff
end

data = [];
err = ~areequal(spc0,spc1);

function w = ord_ones(phi,theta)
w = ones(size(phi));
return;

function w = ord_gaussian(phi,theta)
w = gaussian(theta,0,pi/8);
return;
