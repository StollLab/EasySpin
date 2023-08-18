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

% Isotropic distribution
Exp.Ordering = @ordering_isotropic;
[B,spc1] = pepper(Sys,Exp,Opt);

% Gaussian distribution
Exp.Ordering = @ordering_gaussian;
[B,spc2] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(B,spc0,B,spc1,B,spc2);
  legend('no ordering function','isotropic ordering function','Gaussian');
  legend boxoff
end

ok = areequal(spc0,spc1,1e-10,'rel');

end

function w = ordering_isotropic(beta,gamma)
w = ones(size(beta)).*ones(size(gamma));
end

function w = ordering_gaussian(beta,gamma)
w = gaussian(beta,0,pi/8).*ones(size(gamma));
end
