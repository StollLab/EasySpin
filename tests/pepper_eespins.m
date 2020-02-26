function [ok,data] = test(opt,olddata)

% Regression test: Two coupled electrons
%------------------------------------------------------

% two-electron system, both electrons with orthorhombic g
Sys = struct('S',[1/2 1/2],'g',[2 2.05 2.1; 2.2 2.25 2.3],'HStrain',[1 1 1]*30);
Sys.ee = [1 1 -2]*100; % electron-electron coupling in MHz

% X band conditions
Par = struct('mwFreq',9.5,'Range',[280 350]);
Opt = [];

[x,spc] = pepper(Sys,Par,Opt);

if (opt.Display)
  if isempty(olddata)
    plot(x,spc);
  else
    h = plot(x,spc,'k',x,olddata.spc,'r');
    legend('old','new');
    %set(h(1),'LineWidth',2);
  end
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: two coupled orthorhombic S=1/2');
end

data.spc = spc;

if ~isempty(olddata)
  ok = areequal(spc,olddata.spc,1e-3,'rel');
else
  ok = [];
end
