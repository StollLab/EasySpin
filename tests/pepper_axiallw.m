function [err,data] = test(opt,olddata)

% Regression test: Powder, axial system, isotropic lw

Sys = struct('S',1/2,'g',[2 2 2.2],'lw',1);
Exp = struct('mwFreq',9.5,'Range',[300 350]);
Opt = struct('Verbosity',opt.Verbosity);
[x,spc] = pepper(Sys,Exp,Opt);

if (opt.Display)
  if isempty(olddata)
    plot(x,spc);
  else
    subplot(4,1,[1 2 3]);
    h = plot(x,spc,'k',x,olddata.spc,'r');
    legend('old','new');
    %set(h(1),'LineWidth',2);
    subplot(4,1,4);
    plot(x,spc-olddata.spc);
  end
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Axial system with convolution lw');
end

data.spc = spc;

if ~isempty(olddata)
  err = ~areequal(spc,olddata.spc,1e-4,'rel');
else
  err = [];
end
