function ok = test(opt)

%=======================================================================
% Harmonic setting vs. derivate/integral afterwards
%=======================================================================
% Isotropic convolutional line broadening

Sys = struct('S',1/2,'g',[2 2 2.3]);
Exp = struct('mwFreq',9.7979,'Range',[290 360],'nPoints',10000);
Opt = struct('Verbosity',0,'nKnots',[30 1]);

Sys.lw = 2; % Gaussian broadening

Exp.Harmonic = 0;
[x,y0] = pepper(Sys,Exp,Opt);
y0 = deriv(x,y0);

Exp.Harmonic = 1;
[x,y1] = pepper(Sys,Exp,Opt);

Exp.Harmonic = 2;
[x,y2] = pepper(Sys,Exp,Opt);
y2 = cumtrapz(x,y2);

if (opt.Display)
  subplot(3,1,[1 2]);
  plot(x,y1,x,y2,x,y0);
  legend('Harmonic 1','Harmonic 2, integ','Harmonic 0, diff');
  title('pepper: different harmonics, comparison');
  subplot(3,1,3);
  plot(x,y0-y1,x,y2-y1);
  title('residuals');
  legend('Harmonic 0, diff','Harmonic 2, integ');
end

spread = max(y1)-min(y1);
tolerance = 1e-4;
ok = all(abs(y0-y1)/spread<tolerance) && all(abs(y2-y1)/spread<tolerance);
