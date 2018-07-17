function [err,data] = test(opt,olddata)
% Regression test for MD trajectory-based simulation of EPR spectrum using
% cardamom and ISTOs method (with correlation function approximation)

rng(1)

% Load pre-processed MD frame trajectory
% -------------------------------------------------------------------------

TrajDir = '.\mdfiles\';

load([TrajDir, 'MTSSL_polyAla_traj.mat'])

tScale = 2.5;  % diffusion constant of TIP3P model water molecules in MD 
               % simulations is ~2.5x too high, so we scale the time axis

MD = Traj;
MD.tScale = tScale;
MD.removeGlobal = 0;

MD.DiffGlobal = 6e6;

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------

T = 200e-9;

Sys.Nucs = '14N';

Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10);
Sys.Diff = 3e7;
Sys.lw = [0.1, 0.1];

Par.nTraj = 100;

Par.dt = 2.0e-9;
Par.nSteps = ceil(T/Par.dt);
Par.nOrients = 100;
Par.Model = 'MD';

Exp.mwFreq = 9.4;

Opt.Verbosity = 0;
Opt.FFTWindow = 1;
Opt.Method = 'ISTOs';
Opt.truncate = 30;

[B, spc, ExpectVal,t] = cardamom(Sys,Exp,Par,Opt,MD);
spc = spc/max(spc);

% Plot for comparison
% -------------------------------------------------------------------------

OldDataFile = [TrajDir, 'MTSSL_polyAla_spc_istos.mat'];

if exist(OldDataFile, 'file')>0
  load(OldDataFile)

  diff = abs(spc-spcOld);
  rmsd = sqrt(mean(diff.^2));

  if rmsd < 1e-2
    err = 0;

  else
    err = 1;
    figure

    plot(BOld, spcOld, B, spc)
    ylim([-1.1,1.1])
    ylabel('Im(FFT(M_{+}(t)))')
    xlabel('B (mT)')
    legend('Old','Current')
    hold off
  end
else
  tOld = t;
  ExpectValOld = ExpectVal;
  BOld = B;
  spcOld = spc;
  save(OldDataFile, 'tOld', 'ExpectValOld', 'BOld', 'spcOld')
end

data = [];

end
