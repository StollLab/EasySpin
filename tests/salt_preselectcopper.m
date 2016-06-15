function [err,data] = test(opt,olddata)

% OriPreSelect option, axial g, 63Cu and 1H

Sy = struct('S',1/2,'g',[2 2 2.2],'Nucs','63Cu,1H',...
  'A',[[40 40 400];4+[-1 -1 2]*1.5]);
Sy.lwEndor = 0.02;
Ex = struct('Range',[10 20],'mwFreq',9.5,'ExciteWidth',50);
Ex.Field = Ex.mwFreq*1e9*planck/bmagn/2.05*1e3;
Op = struct('Verbosity',opt.Verbosity,'nKnots',91);

Op.OriPreSelect = 0; [x,y0] = salt(Sy,Ex,Op);
Op.OriPreSelect = 1; [x,y1] = salt(Sy,Ex,Op);

if (opt.Display)
  plot(x,y0,'b',x,y1,'r');
end

RelativeError = abs(y0(:)-y1(:))/max(y0(:));
err = any(RelativeError>0.1);
data = [];
