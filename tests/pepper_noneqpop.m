function [err,data] = test(opt,olddata)

% Regression test: Non-equilibrium populations
%===============================================

% Spin system and experimental parameters
Sys = struct('S',1,'g',2,'lw',0.3,'D',300);
Exp = struct('mwFreq',9.5,'Range',[325 355],'Harmonic',0);

% User-specified population vector
Sys.Pop = [0.85 1 0.95];

% Simulation options
Opt = struct;

[x,spc] = pepper(Sys,Exp,Opt);


if opt.Display
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
  title('pepper: Nonequilibrum populations');
end

data.spc = spc;

if ~isempty(olddata)
  err = ~areequal(spc,olddata.spc,1e-4,'rel');
else
  err = [];
end
