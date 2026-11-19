function ok = test() %#ok<FNDEF>

% With Opt.IdealPulses, finite pulses are treated as ideal, with a
% Gaussian excitation window of FWHM ExciteWidth = 1/tp. The excitation
% weight is the fraction of the inhomogeneous EPR line (FWHM from lwAnisoG)
% within the window, ExciteWidth/sqrt(ExciteWidth^2+lw^2). Checked for the
% S=1/2 path and the general path (a second nucleus in the core system,
% via Opt.Nuclei; S > 1/2 does not support Opt.IdealPulses).

Sys{1} = struct('S',1/2,'g',2,'Nucs','1H','A',0.01);
Sys{2} = struct('S',1/2,'g',2,'Nucs','1H,1H','A',[0.01; 0.01]);
Nuclei = {[],1};
Exp = struct('Sequence','2pESEEM','Field',350,'tau',0.1,'dt',0.01,'nPoints',4);
Exp.SampleFrame = [0 0 0];
ExciteWidth = 20; % MHz
Exp.tp = 1/ExciteWidth;  % µs
lwInh = 30; % MHz

nu0 = 2*bmagn*Exp.Field*1e-3/planck/1e9; % EPR frequency, GHz

tol = 1e-4;
for k = 1:2
  I0 = echo(Sys{k},Exp,nu0,0,0,Nuclei{k});
  Ihalf = echo(Sys{k},Exp,nu0,ExciteWidth/2,0,Nuclei{k});
  Iinh = echo(Sys{k},Exp,nu0,0,lwInh,Nuclei{k});
  ok(2*k-1) = areequal(Ihalf/I0,0.5,tol,'rel');
  ok(2*k) = areequal(Iinh/I0,ExciteWidth/sqrt(ExciteWidth^2+lwInh^2),tol,'rel');
end

end

% Echo amplitude for given mw offset (MHz) and inhomogeneous width (MHz),
% with the nuclei of the nuclear subspace (empty: all)
function I = echo(Sys,Exp,nu0,offset,lw,Nuclei)
Sys.lwAnisoG = lw;
Exp.mwFreq = nu0 + offset/1e3;
Opt.IdealPulses = true;
Opt.Method = 'subspace';
Opt.Nuclei = Nuclei;
[~,y] = saffron(Sys,Exp,Opt);
I = real(y(1));
end
