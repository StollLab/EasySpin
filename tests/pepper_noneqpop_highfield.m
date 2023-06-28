function [ok,data] = test(opt,olddata)

% Regression test: Non-equilibrium populations
%                  (high-field)
%===============================================

% Spin system and experimental parameters
Sys = struct('S',1,'g',2,'lw',0.3,'D',[300 30]);
Exp = struct('mwFreq',9.5,'Range',[325 355],'Harmonic',0);

% Simulation options
Opt = struct;

% User-specified high-field population vector
Sys.initState = {[0 1 0],'eigen'};

[x,spc] = pepper(Sys,Exp,Opt);

% Short-cut for T0-populated triplet state
Sys.initState = 'T0';
[~,spc_] = pepper(Sys,Exp,Opt);

% Frequency sweep
Sys.lw = unitconvert(Sys.lw,'mT->MHz');
Expf = struct('Field',340);

[f,spcf] = pepper(Sys,Expf,Opt);

% Short-cut for T0-populated triplet state
Sys.initState = 'T0';
[~,spcf_] = pepper(Sys,Expf,Opt);

if opt.Display
  if isempty(olddata)
    subplot(2,1,1)
    plot(x,spc);
    xlabel('magnetic field (mT)')
    subplot(2,1,2)
    plot(f,spcf);
    xlabel('microwave frequency (GHz)')
  else
    subplot(4,2,[1 3 5]);
    plot(x,spc,'k',x,olddata.spc,'r');
    legend('new','old');
    subplot(4,2,7);
    plot(x,spc-olddata.spc);
    xlabel('magnetic field (mT)');
    ylabel('intensity (a.u.)');
    subplot(4,2,[2 4 6]);
    plot(f,spcf,'k',f,olddata.spcf,'r');
    legend('new','old');
    subplot(4,2,8);
    plot(f,spcf-olddata.spcf);
    xlabel('microwave frequency (GHz)');
    ylabel('intensity (a.u.)');
  end
  sgtitle('pepper: Nonequilibrum populations recombination triplet','FontSize',9,'FontWeight','bold');
end

data.spc = spc;
data.spcf = spcf;

ok = areequal(spc,spc_,1e-12,'abs') && areequal(spcf,spcf_,1e-12,'abs');
if ~isempty(olddata)
  ok = [ok areequal(spc,olddata.spc,1e-4,'rel') && areequal(spcf,olddata.spcf,1e-4,'rel')];
end
