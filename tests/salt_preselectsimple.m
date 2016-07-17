function [err,data] = test(opt,olddata)

% OriPreSelect option, axial g and 1H

Sy = struct('S',1/2,'g',[2 2 2.2],'Nucs','1H','A',4+[-1 -1 2]*1.5);
Sy.lwEndor = 0.02;
Ex = struct('Range',[8 18],'mwFreq',9.5,'ExciteWidth',100);
g = linspace(min(Sy.g),max(Sy.g),10);
Fields = Ex.mwFreq*1e9*planck/bmagn./g*1e3;
Op = struct('Verbosity',opt.Verbosity,'nKnots',61);

for iField = 1:numel(Fields)
  Ex.Field = Fields(iField);
  Op.OriPreSelect = 0;
  [x,y0(iField,:)] = salt(Sy,Ex,Op);
  Op.OriPreSelect = 1;
  [x,y1(iField,:)] = salt(Sy,Ex,Op);
end

Residual = y1-y0;

if (opt.Display)
  subplot(3,1,1); plot(x,y0,'k'); legend('without OriPreSelect');
  title('OriPreSelect option, simple test');
  subplot(3,1,2); plot(x,y1,'g'); legend('with OriPreSelect');
  xlabel('frequency [MHz]');
  subplot(3,1,3); plot(x,Residual,'r'); legend('Residual');
end

RelativeError = max(abs(Residual(:)))/max(y0(:));
err = any(RelativeError>0.02);

data = [];
