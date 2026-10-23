function ok = test()

% Exp.ExciteWidth is the FWHM of the Gaussian excitation window, and the
% excitation weight is the fraction of the inhomogeneous EPR line (FWHM
% from lwAnisoG) within the window, ExciteWidth/sqrt(ExciteWidth^2+lw^2).

Sys = struct('S',1/2,'g',2,'Nucs','1H','A',1e-3); % small A: both EPR lines at the same offset
Exp = struct('Field',350,'SampleFrame',[0 0 0]);
Exp.ExciteWidth = 20; % MHz
lwInh = 30; % MHz

nu0 = 2*bmagn*Exp.Field*1e-3/planck/1e9; % EPR frequency, GHz

% ENDOR intensity for given mw offset (MHz) and inhomogeneous width
intensity = @(fun,offset,lw) firstint(fun,setfield(Sys,'lwAnisoG',lw),...
  setfield(Exp,'mwFreq',nu0+offset/1e3));

tol = 1e-5;
funs = {@endorfrq,@endorfrq_perturb};
for k = 1:numel(funs)
  f = funs{k};
  I0 = intensity(f,0,0);
  Ihalf = intensity(f,Exp.ExciteWidth/2,0);
  Iinh = intensity(f,0,lwInh);
  ok(2*k-1) = areequal(Ihalf/I0,0.5,tol,'rel');
  ok(2*k) = areequal(Iinh/I0,Exp.ExciteWidth/sqrt(Exp.ExciteWidth^2+lwInh^2),tol,'rel');
end

end

function I = firstint(fun,Sys,Exp)
[~,Int] = fun(Sys,Exp);
I = Int(1);
end
