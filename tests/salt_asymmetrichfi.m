function [err,data] = test(opt,olddata)

% Tests if asymmetric A tensor work

Sy = struct('S',1/2,'Nucs','1H','g',[2.2 2.1 2],...
  'A',2*[-1.2,-0.8,2],'lwEndor',[0,0.01]);
Ex = struct('Range',[10 18],'Field',326.5,...
  'nPoints',2^12,'mwFreq',9.7);
Si = struct('Threshold',1e-5,'nKnots',[20 5],'Intensity','on',...
  'Enhancement','off','Verbosity',opt.Verbosity);

Sy.A = diag(Sy.A);
[a,b1] = salt(Sy,Ex,Si);
Sy.A(1,2) = 5;
Sy.A(2,1) = -5;
[a,b2] = salt(Sy,Ex,Si);

if opt.Display
  subplot(3,1,[1 2]);
  plot(a,b1,'b',a,b2,'r');
  title('Tilted tensors, Ci symmetry');
  legend('symmetric A','asymmetric A');
  xlabel('frequency [MHz]');
  if ~isempty(olddata)
    subplot(3,1,3);
    plot(a,b1-olddata.b1,'b',a,b2-olddata.b2,'r');
  end
end

data.b1 = b1;
data.b2 = b2;

if isempty(olddata)
  err = [];
else
  ok = areequal(olddata.b1,b1,1e-10,'rel') && areequal(olddata.b2,b2,1e-10,'rel');
  err = ~ok;
end
