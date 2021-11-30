function [ok,data] = test(opt,olddata)

% Regression test: Spin-correlated radical pair
%===============================================

% Spin system
Sys = struct('S',[1/2 1/2],'g',[2.00 1.99],'HStrain',1);
Sys.ee = [1 1 -2]*5; % electron-electron coupling in MHz

% Experimental parameters
Exp = struct('mwFreq',9.5,'Range',[337 343],'Harmonic',0);

% User-specified population vector
S = 1/sqrt(2)*[0 1 -1 0];
Sys.initState = S'*S;

% Simulation options
Opt = struct;

[x,spc] = pepper(Sys,Exp,Opt);

% With hyperfine coupling
SysHF = Sys;
SysHF.Nucs = '1H';
SysHF.A = [2 2 2 0 0 0]; % MHz

[x,spcHF] = pepper(SysHF,Exp,Opt);

if opt.Display
  if isempty(olddata)
    plot(x,spc,x,spcHF);
  else
    subplot(4,2,[1 3 5]);
    h1 = plot(x,spc,'k',x,olddata.spc,'r');
    legend('old','new');
    %set(h1(1),'LineWidth',2);
    subplot(4,2,7);
    plot(x,spc-olddata.spc);
    subplot(4,2,[2 4 6]);
    h2 = plot(x,spcHF,'k',x,olddata.spcHF,'r');
    legend('old','new');
    %set(h2(1),'LineWidth',2);
    subplot(4,2,8);
    plot(x,spcHF-olddata.spcHF);
  end
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Spin-correlated radical pair');
end

data.spc = spc;
data.spcHF = spcHF;

if ~isempty(olddata)
  ok = areequal(spc,olddata.spc,1e-4,'rel') && areequal(spcHF,olddata.spcHF,1e-4,'rel');
else
  ok = [];
end
