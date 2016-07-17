function [err,data] = test(opt,olddata)

% Rotational invariance

Sys.Nucs = '1H';
Sys.g = [2 2.05 2.1];
Sys.A = 2*[-1.2,-0.8,2];
Sys.lwEndor = 0.1;

Exp.Field = 326.5;
Exp.Range = [10 20];

Opt.Threshold = 1e-5;
Opt.nKnots = [30 5];
Opt.Verbosity = opt.Verbosity;

[rf,b0] = salt(Sys,Exp,Opt);
Sys.AFrame = pi/180*[-34 19 -45];
[rf,b1] = salt(Sys,Exp,Opt);
Sys.AFrame = pi/180*[55 -72 -5];
[rf,b2] = salt(Sys,Exp,Opt);

if (opt.Display)
  subplot(3,1,[1 2]);
  plot(rf,b0,'b',rf,b1,'r',rf,b2,'g');
  title('Rotation invariance, rhombic A');
  legend('rhombic A along xyz (D2h)','A tilted (Ci)','A tilted again (Ci)');
  legend boxoff
  xlabel('frequency [MHz]');
  if ~isempty(olddata)
    subplot(3,1,3);
    plot(rf,b0-olddata.b0,'b',rf,b1-olddata.b1,'r',rf,b2-olddata.b2,'g');
    title('Regresion test, residuals');
  end
end

data.b0 = b0;
data.b1 = b1;
data.b2 = b2;

if isempty(olddata)
  err = [];
else
  thr = max(b0)*1e-4;
  ok = areequal(olddata.b0,b0,thr) & ...
       areequal(olddata.b1,b1,thr) & ...
       areequal(olddata.b2,b2,thr);
  err = ~ok;
end
