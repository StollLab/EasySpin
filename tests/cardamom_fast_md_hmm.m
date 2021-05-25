function [ok,data] = test(opt,olddata)
% Regression test for MD trajectory-based simulation of EPR spectrum using
% cardamom and MSM method for nitroxides

rng(1)

% Load pre-processed MD frame trajectory
% -------------------------------------------------------------------------

load('.\mdfiles\MTSSL_polyAla_traj.mat')
MD = Traj;

% Correct array sizes such that nTraj is last dimension
f = fieldnames(MD);
for k = 1:numel(f)
  fn = f{k};
  if ndims(MD.(fn))~=4, continue; end
  MD.(fn) = permute(MD.(fn),[1 2 4 3]);
end

tScale = 2.5;  % diffusion constant of TIP3P model water molecules in MD 
               % simulations is ~2.5x too high, so we scale the time axis

MD.dt = MD.dt*tScale;
MD.removeGlobal = false;
MD.tLag = 100e-12*tScale;
MD.nStates = 20;

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------
Sys.g = [2.009, 2.006, 2.002];
Sys.Nucs = '14N';
Sys.A = mt2mhz([6, 36]/10);
Sys.lw = [0.1, 0.1];

Par.dtSpatial = MD.tLag;
T = 250e-9;
Par.dtSpin = Par.dtSpatial;
Par.nSteps = ceil(T/Par.dtSpin);

Par.Model = 'MD-HMM';
Par.nTraj = 100;
Par.nOrients = 100;

Exp.mwFreq = 9.4;

Opt.Verbosity = 0;
Opt.FFTWindow = 1;
Opt.Method = 'fast';
Opt.nTrials = 2;

[B, spc] = cardamom(Sys,Exp,Par,Opt,MD);
spc = spc/max(spc);

data.spc = spc;

if ~isempty(olddata)
  ok = areequal(olddata.spc,spc,1e-2,'abs');
else
  ok = [];
end

if opt.Display
  if ~isempty(olddata)
    plot(B,olddata.spc,B,data.spc);
    legend('old','new');
    legend boxoff
    title(mfilename,'Interpreter','none');
    axis tight
  end
end

end
