function ok = test()

% Exp.ExciteWidth is the FWHM of the Gaussian excitation window, and the
% excitation weight is the fraction of the inhomogeneous EPR line (FWHM
% from lwAnisoG) within the window, ExciteWidth/sqrt(ExciteWidth^2+lw^2).
% Checked for the S=1/2 path and the general path (S=1).

Sys{1} = struct('S',1/2,'g',2,'Nucs','1H','A',0.01);
Sys{2} = struct('S',1,'g',2,'Nucs','1H','A',0.01);  % S=1, D=0: both EPR transitions at nu0
Exp = struct('Sequence','2pESEEM','Field',350,'tau',0.1,'dt',0.01,'nPoints',4);
Exp.SampleFrame = [0 0 0];
Exp.ExciteWidth = 20; % MHz
lwInh = 30; % MHz

nu0 = 2*bmagn*Exp.Field*1e-3/planck/1e9; % EPR frequency, GHz

tol = 1e-4;
for k = 1:2
  I0 = echo(Sys{k},Exp,nu0,0,0);
  Ihalf = echo(Sys{k},Exp,nu0,Exp.ExciteWidth/2,0);
  Iinh = echo(Sys{k},Exp,nu0,0,lwInh);
  ok(2*k-1) = areequal(Ihalf/I0,0.5,tol,'rel');
  ok(2*k) = areequal(Iinh/I0,Exp.ExciteWidth/sqrt(Exp.ExciteWidth^2+lwInh^2),tol,'rel');
end

end

% Echo amplitude for given mw offset (MHz) and inhomogeneous width (MHz)
function I = echo(Sys,Exp,nu0,offset,lw)
Sys.lwAnisoG = lw;
Exp.mwFreq = nu0 + offset/1e3;
[~,y] = saffron(Sys,Exp);
I = real(y(1));
end
