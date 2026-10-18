function ok = test(opt)

% Check that the Lorentzian line width of a frequency sweep is in MHz for all
% accumulation methods

Sys.g = 2;
Sys.lw = [0 5];  % MHz
Exp.Field = 340;
Exp.mwRange = [9.45 9.6];
Exp.nPoints = 2000;
Exp.Harmonic = 0;

Opt.AccumMethod = 'binning';
[x,spc0] = garlic(Sys,Exp,Opt);
Opt.AccumMethod = 'explicit';
[x,spc1] = garlic(Sys,Exp,Opt);
Opt.AccumMethod = 'template';
[x,spc2] = garlic(Sys,Exp,Opt);

% Exact Lorentzian, normalized to its maximum
nu0 = bmagn*2*Exp.Field*1e-3/planck/1e9;  % GHz
spcRef = lorentzian(x,nu0,Sys.lw(2)/1e3);
spcRef = spcRef/max(spcRef);

ok(1) = areequal(spc0/max(spc0),spcRef,2e-2,'abs');
ok(2) = areequal(spc1/max(spc1),spcRef,2e-2,'abs');
ok(3) = areequal(spc2/max(spc2),spcRef,2e-2,'abs');

if opt.Display
  plot(x,spcRef,x,spc0/max(spc0),'--',x,spc1/max(spc1),':',x,spc2/max(spc2),'-.');
  legend('exact','binning','explicit','template');
end

end
