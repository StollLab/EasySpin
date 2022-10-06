function [ok,data] = test(opt,olddata)

% Regression test: ISC triplet state with axial ZFS
%===================================================

% Spin system and experimental parameters
Sys = struct('S',1,'g',2,'lw',0.3,'D',500);
Exp = struct('mwFreq',9.5,'Range',[310 370],'Harmonic',0);

% User-specified population vector
pop = [0.3 0.6 0.1]; % populations of zero-field levels

% Get zero-field states and order in terms of energy
[F,~,~,~] = ham(Sys);
[ZFStates,ZFEnergies] = eig(F);
[ZFEnergies,idx] = sort(real(diag(ZFEnergies)));
ZFStates = ZFStates(:,idx);

% Correct zero-field states for S=1 and axial D
if ZFEnergies(2)==ZFEnergies(3)
  % Manual zero-field states (D>0)
  v1 = ZFStates(:,2);
  v2 = ZFStates(:,3);
  ZFStates(:,2) = (v1-v2)/sqrt(2);
  ZFStates(:,3) = (v1+v2)/sqrt(2);
elseif ZFEnergies(2)==ZFEnergies(1)
  % Manual zero-field states (D<0)
  v1 = ZFStates(:,1);
  v2 = ZFStates(:,2);
  ZFStates(:,2) = (v1-v2)/sqrt(2);
  ZFStates(:,1) = (v1+v2)/sqrt(2);
end

% Convert population vector to density matrix
Sys.initState = ZFStates*diag(pop)*ZFStates';

% Simulation options
Opt = struct;

[x,spc] = pepper(Sys,Exp,Opt);

% Frequency sweep
Sys.lw = mt2mhz(Sys.lw);
Expf = struct('Field',340,'mwRange',[8.5 10.5]);

[f,spcf] = pepper(Sys,Expf,Opt);

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
    legend('old','new');
    subplot(4,2,7);
    plot(x,spc-olddata.spc);
    xlabel('magnetic field (mT)');
    ylabel('intensity (a.u.)');
    subplot(4,2,[2 4 6]);
    plot(f,spcf,'k',f,olddata.spcf,'r');
    legend('old','new');
    subplot(4,2,8);
    plot(f,spcf-olddata.spcf);
    xlabel('microwave frequency (GHz)');
    ylabel('intensity (a.u.)');
  end
  sgtitle('pepper: Nonequilibrum populations ISC triplet','FontSize',9,'FontWeight','bold');
end

data.spc = spc;
data.spcf = spcf;

if ~isempty(olddata)
  ok = areequal(spc,olddata.spc,1e-4,'rel') && areequal(spcf,olddata.spcf,1e-4,'rel');
else
  ok = [];
end
