function [ok,data] = test(opt,olddata)

% Triplet pair or quintet state formed by 
% singlet fission
% (direct test for comparison of Sys.S = [1 1]
%  and Sys.S = 2 and regression test)
%===============================================

% Spin system and experimental parameters
D = 1400; % MHz
SysTT = struct('S',[1 1],'lw',2,'D',[D; D],'J',1e7);
SysQ = struct('S',2,'lw',2,'D',D/3);
Exp = struct('mwFreq',9.5,'Range',[300 380],'Harmonic',0);

% Simulation options
Opt = struct;

% User-specified high-field population vectors
SysTT.initState = {[0 0 0 0 0 0 1 0 0],'eigen'};
SysQ.initState = {[0 0 1 0 0],'eigen'};

[x,spcTT] = pepper(SysTT,Exp,Opt);
[~,spcQ] = pepper(SysQ,Exp,Opt);

if opt.Display
  if isempty(olddata)
    plot(x,spcTT,x,spcQ);
    xlabel('magnetic field (mT)')
  else
    subplot(4,2,[1 3 5]);
    plot(x,spcTT,'k',x,olddata.spcTT,'r');
    legend('new','old');
    subplot(4,2,7);
    plot(x,spcTT-olddata.spcTT);
    xlabel('magnetic field (mT)');
    ylabel('intensity (a.u.)');
    subplot(4,2,[2 4 6]);
    plot(x,spcQ,'k',x,olddata.spcQ,'r');
    legend('new','old');
    subplot(4,2,8);
    plot(x,spcQ-olddata.spcQ);
    xlabel('magnetic field (mT)');
    ylabel('intensity (a.u.)');
  end
  sgtitle('pepper: triplet pair/quintet state from singlet fission','FontSize',9,'FontWeight','bold');
end

data.spcTT = spcTT;
data.spcQ = spcQ;

ok = areequal(spcTT,spcQ,1e-3,'rel');
if ~isempty(olddata)
  ok = [ok areequal(spcTT,olddata.spcTT,1e-4,'rel') && areequal(spcQ,olddata.spcQ,1e-4,'rel')];
end
