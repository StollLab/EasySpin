function [ok,data] = test(opt,olddata)

%=======================================================================
% Powder, orthorhombic system, isotropic lw
%=======================================================================
Sys = struct('S',1/2,'g',[2 2.1 2.2],'lw',1);
Exp = struct('mwFreq',9.5,'Range',[300 350]);
Opt = struct();
[x,y1] = pepper(Sys,Exp,Opt);

if (opt.Display)
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(x,y1,x,olddata.y1);
    ylabel('intensity (arb.u.)');
    title('Orthorhombic system with convolution lw');
    subplot(3,1,3);
    plot(x,y1-olddata.y1,'r');
    xlabel('magnetic field [mT]');
    title('Residuals');
  end
end

data.y1 = y1;

if ~isempty(olddata)
  ok = areequal(y1,olddata.y1,1e-3,'rel');
else
  ok = [];
end

